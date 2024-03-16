import re
import pysam
import pandas as pd
import numpy as np
import subprocess
import shlex
import gzip

rng = np.random.default_rng(1)

# R1 is G to A mutated in unmethylated C
R1_CUT_SITES = [
    'CATG',  # NlaIII
    'CATA',  # NlaIII
    'GATC',  # DpnII or MboI
    'AATC'  # DpnII or MboI
]

r1_cut = re.compile("(" + "|".join(R1_CUT_SITES) + ")")

# R2 is C to T mutated in unmethylated C
R2_CUT_SITES = [
    'CATG',  # NlaIII
    'TATG',  # NlaIII
    'GATC',  # DpnII or MboI
    'GATT'  # DpnII or MboI
]

r2_cut = re.compile("(" + "|".join(R2_CUT_SITES) + ")")

ref_cut = re.compile("(GATC|CATG)")

def illegal_overlap(seg_keys):
    illegal_overlap = False
    # If chimera, determine if sequential overlaps condition is violated
    if len(seg_keys) > 2:
        for i in range(len(seg_keys) - 1):
            key = seg_keys[i]
            next_key = seg_keys[i+1]
            if i == 0:
                previous_key = (-1, -1, 0) # Dummy key
            else:
                previous_key = seg_keys[i-1]
        
            # Current key is contained within previous key
            if key[0] >= previous_key[0] and key[1] <= previous_key[1]:
                illegal_overlap = True
                break
            # Current key contains previous key
            if key[0] <= previous_key[0] and key[1] >= previous_key[1]:
                illegal_overlap = True
                break
            # Current key is contained within next key
            elif key[0] >= next_key[0] and key[1] <= next_key[1]:
                illegal_overlap = True
                break
            # Current key contains next key
            elif key[0] <= next_key[0] and key[1] >= next_key[1]:
                illegal_overlap = True
                break
        
            # Check for overlaps of current key with all other keys that are 
            # not the previous or next key
            for j in range(len(seg_keys)):
                test_seg_key = seg_keys[j]
                if test_seg_key == previous_key:
                    continue
                elif test_seg_key == next_key:
                    continue
                elif test_seg_key == key:
                    continue
                elif test_seg_key[0] < key[0] < test_seg_key[1]:
                    illegal_overlap = True
                    break
                elif test_seg_key[0] < key[1] < test_seg_key[1]:
                    illegal_overlap = True
                    break

    elif len(seg_keys) == 2:
        key = seg_keys[0]
        next_key = seg_keys[1]

        # Key is contained within next key
        if key[0] >= next_key[0] and key[1] <= next_key[1]:
            illegal_overlap = True
        # Key contains next key
        elif key[0] <= next_key[0] and key[1] >= next_key[1]:
            illegal_overlap = True

    return illegal_overlap

def contact_filter(algn1, 
                   algn2, 
                   pair_index,
                   R1_cs_keys,
                   R2_cs_keys
                  ):
    if not algn1["is_mapped"] or not algn1["is_unique"]:
        return "na"
    if not algn2["is_mapped"] or not algn2["is_unique"]:
        return "na"
    contact_type = pair_index[1]
    
    if contact_type == "R1": # Pairtools reports 5' fragment before 3' fragment
        idx5 = algn1["idx"]
        idx3 = algn2["idx"]
        cs_key = (idx5, idx3)
        if cs_key not in R1_cs_keys:
            return "chimera"
        if not R1_cs_keys[cs_key]:
            return "chimera"
    elif contact_type == "R2": # Pairtools reports 3' fragment before 5' fragment
        idx5 = algn2["idx"]
        idx3 = algn1["idx"]
        cs_key = (idx5, idx3)
        if cs_key not in R2_cs_keys:
            return "chimera"
        if not R2_cs_keys[cs_key]:
            return "chimera"
    return "contact"

def get_methylation_tags(pairs, mate, is_forward):

    ZC = 0
    ZR = 0
    NM = 0
    possible_conv = 0
    
    for pair in pairs:
        #print(pair)
        if pair["cigar"] == 0:
            if pair["ref_nuc"].islower():
                NM += 1
                if pair["ref_nuc"] == "c" and pair["read_nuc"] == "T":
                    possible_conv += 1
                elif pair["ref_nuc"] == "g" and pair["read_nuc"] == "A":
                    possible_conv += 1
            if mate == "2":
                if is_forward:
                    if pair["ref_nuc"] == "c" and pair["read_nuc"] == "T":
                        ZC += 1
                    elif pair["ref_nuc"] == "C" and pair["read_nuc"] == "C":
                        ZR += 1
                else:
                    if pair["ref_nuc"] == "g" and pair["read_nuc"] == "A":
                        ZC += 1
                    elif pair["ref_nuc"] == "G" and pair["read_nuc"] == "G":
                        ZR += 1
            elif mate == "1":
                if is_forward:
                    if pair["ref_nuc"] == "g" and pair["read_nuc"] == "A":
                        ZC += 1
                    elif pair["ref_nuc"] == "G" and pair["read_nuc"] == "G":
                        ZR += 1
                else:
                    if pair["ref_nuc"] == "c" and pair["read_nuc"] == "T":
                        ZC += 1
                    elif pair["ref_nuc"] == "C" and pair["read_nuc"] == "C":
                        ZR += 1
        elif pair["cigar"] in [1, 2]:
            NM += 1
        #print(ZC, ZR, NM)
    
    NM = NM - possible_conv

    #print(possible_conv)

    return ZC, ZR, NM

def get_md_tag(pairs):
    MD = ""
    match_count = 0
    deletion = False

    for i in range(len(pairs)):
        pair = pairs[i]
        cigar_id = pair["cigar"]
        if cigar_id == 0:
            deletion = False
            if pair["ref_nuc"].isupper():
                match_count += 1
            else:
                MD += str(match_count)
                MD += pair["ref_nuc"].upper()
                match_count = 0
        elif cigar_id == 2:
            if deletion == False:
                deletion = True
                MD += str(match_count)
                MD += "^"
                match_count = 0
            MD += pair["ref_nuc"].upper()
    MD += str(match_count)

    return MD

def softclip_pairs(pairs):
    new_pairs = []
    for i in range(len(pairs)):
        if pairs[i]["cigar"] in [0,1]:
            new_algn = {
                "cigar" : 4,
                "read_pos" : None,
                "ref_pos" : None,
                "ref_nuc" : None,
                "read_nuc" : None,
                "read_qual" : None,
                #"original_sequence_nuc" : None
            }
            new_pairs.append(new_algn)
        elif pairs[i]["cigar"] in [2]:
            continue
        else:
            new_pairs.append(pairs[i])
        
    return new_pairs            

def create_trimmed_mate(read, original_span, adjusted_span, pairs, mate, header):

    if original_span == adjusted_span:
        return read, pairs
        
    pos_5, index_5, read_5 = get_5_cut_point(pairs, adjusted_span[0], read.is_forward)
    # Index is at end of range; subtract 1 to get actual final index
    pos_3, index_3, read_3 = get_3_cut_point(pairs, adjusted_span[1] - 1, read.is_forward)

    if read.is_forward:
        start_index = index_5
        end_index = index_3
        pos = pos_5
    else:
        start_index = index_3
        end_index = index_5
        pos = pos_3

    if start_index >= end_index:
        return None, None
    
    new_pairs = pairs[start_index:end_index]
    if not read.is_secondary:
        left_pairs = softclip_pairs(pairs[:start_index])
        right_pairs = softclip_pairs(pairs[end_index:])
        new_pairs = left_pairs + new_pairs + right_pairs

    if read.is_secondary:
        read_seq = ""
        read_qual = ""
        for pair in new_pairs:
            if pair["cigar"] in [0, 1]:
                read_seq += pair["read_nuc"]
                read_qual += pair["read_qual"]
    else:
        read_seq = read.query_sequence
        read_qual = pysam.qualities_to_qualitystring(read.query_qualities)

    cigartuples = build_cigar(read, pairs, start_index, end_index)
    new_ZC, new_ZR, new_NM = get_methylation_tags(new_pairs, mate, read.is_forward)
    new_MD = get_md_tag(new_pairs)

    new_tags = {
        "ZC" : new_ZC,
        "ZR" : new_ZR,
        "NM" : new_NM,
        "MD" : new_MD
    }

    tags = read.get_tags()

    for i in range(len(tags)):
        tag = tags[i]
        if tag[0] in new_tags:
            tags[i] = (tag[0], new_tags[tag[0]])

    real_trim_5 = read_5 - original_span[0]
    real_trim_3 = original_span[1] - read_3 + 1 

    tags.append(("ZU", real_trim_5))
    tags.append(("ZD", real_trim_3))
    
    new_read = pysam.AlignedSegment(header=pysam.AlignmentHeader().from_dict(header))
    new_read.query_name = read.query_name
    new_read.query_sequence = read_seq
    new_read.flag = read.flag
    new_read.reference_id = read.reference_id
    new_read.reference_start = pos
    new_read.mapping_quality = read.mapping_quality
    new_read.cigar = cigartuples
    new_read.query_qualities = pysam.qualitystring_to_array(read_qual)
    new_read.tags = tags

    return new_read, new_pairs

def get_loc(read):
    cigar = read.cigartuples
    if not cigar:
        return (0, 0, read.mapping_quality)
    #clip_5 = 0
    #clip_3 = 0
    if read.is_forward:
        clip_5 = cigar[0][1] if cigar[0][0] in [4, 5] else 0
        clip_3 = cigar[-1][1] if cigar[-1][0] in [4, 5] else 0
    else:
        clip_5 = cigar[-1][1] if cigar[-1][0] in [4, 5] else 0
        clip_3 = cigar[0][1] if cigar[0][0] in [4, 5] else 0

    start = clip_5
    end = read.get_tag("XL") - clip_3

    return (start, end, read.mapping_quality)

def reverse_complement(seq):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return "".join(nn[n] for n in reversed(seq))

def complement(nuc):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return nn[nuc.upper()]

def alignment_info(read, seg, original_sequence):
    pairs = read.get_aligned_pairs(with_seq = True)
    cigar = read.cigartuples
    read_seq = read.query_sequence
    read_qualities = pysam.qualities_to_qualitystring(read.query_qualities)
    
    if read.is_forward:
        read_start = seg[0]
        increment = 1
    else:
        read_start = seg[1] - 1 # Due to indexing
        increment = -1
    
    full_cigar = []
    for i in cigar:
        if i[0] != 5: # Do not include hard clipping (on secondary alignments)
            full_cigar += [i[0]] * i[1]

    new_pairs = []

    seq_idx = 0
    for i in range(len(pairs)):
        pair = pairs[i]
        cigar_val = full_cigar[i]
        algn_pos = pair[0]
        ref_pos = pair[1]
        ref_nuc = pair[2]
        read_nuc = None
        qual_val = None
        read_pos = None
        original_read_nuc = None
        
        if cigar_val in [0,1, 4]:
            read_nuc = read_seq[seq_idx]
            qual_val = read_qualities[seq_idx]
            if cigar_val != 4:
                read_pos = read_start 
                #original_read_nuc = original_sequence[read_pos]
                read_start += increment
            seq_idx += 1            

        pair_dict = {
            "cigar":cigar_val,
            "algn_pos": algn_pos,
            "read_pos": read_pos,
            "ref_pos": ref_pos,
            "ref_nuc": ref_nuc,
            "read_nuc" : read_nuc,
            "read_qual" : qual_val
        }
        new_pairs.append(pair_dict)

    return new_pairs

def overlap_reference_alignment(aligned_pairs, overlap_start, overlap_end, is_forward):
    seq = ""
    collect = False
    if is_forward:
        iterator = range(len(aligned_pairs))
    else:
        iterator = reversed(range(len(aligned_pairs)))

    for i in iterator:
        pair = aligned_pairs[i]
        if pair["read_pos"] == overlap_start:
            collect = True
        if collect:
            nuc = pair["ref_nuc"]
            if pair['cigar'] == 1:
                nuc = "-" # Add gap for insertion
            elif not is_forward:
                nuc = complement(nuc)
            seq += nuc.upper()
            
        if pair["read_pos"] == overlap_end-1:
            break

    return seq

def get_cut_site_spans(site, seq):
    spans = [i.span() for i in re.finditer(site, seq)]
    return spans
    
def adjust_split(pairs5, seg5, read5,
                 pairs3, seg3, read3,
                 original_sequence, mate):
    
    overlap = seg5[1] - seg3[0]
    cut_site_present = False
    
    if read5.mapping_quality > read3.mapping_quality:
        r5_cover = True
    elif read5.mapping_quality == read3.mapping_quality:
        r5_len = len(read5.query_alignment_sequence)
        r3_len = len(read3.query_alignment_sequence)
        if r5_len > r3_len:
            r5_cover = True
        elif r5_len == r3_len:
            r5_cover = rng.choice([True, False])
        else:
            r5_cover = False
    else:
        r5_cover = False

    if 4 <= overlap  < 8 : # Valid overlap between alignments
        r5_ref_overlap_seq = overlap_reference_alignment(pairs5, seg3[0], seg5[1], read5.is_forward)
        r3_ref_overlap_seq = overlap_reference_alignment(pairs3, seg3[0], seg5[1], read3.is_forward)

        read_overlap_seq = original_sequence[seg3[0]:seg5[1]]

        # Handle insertions/deletions
        if len(r5_ref_overlap_seq) == len(r3_ref_overlap_seq) and \
            "-" not in r5_ref_overlap_seq and "-" not in r3_ref_overlap_seq:
            if mate == "1":
                read_cuts = get_cut_site_spans(r1_cut, read_overlap_seq)
            elif mate == "2":
                read_cuts = get_cut_site_spans(r2_cut, read_overlap_seq)
                
            r5_cut = get_cut_site_spans(ref_cut, r5_ref_overlap_seq)
            r3_cut = get_cut_site_spans(ref_cut, r3_ref_overlap_seq)

            # Given the multiple possibilities of cut sites that can be empirically observed,
            # at least one must be observed. Only one can be observed in r5 and r3 reference,
            # and it must be the same for each.
            if len(read_cuts) > 0 and (r5_cut == r3_cut) and (len(r5_cut) == len(r3_cut) == 1):
                
                cut_span = r5_cut[0]
                if r5_cover:
                    trim_site = seg3[0] + cut_span[1]
                else:
                    trim_site = seg3[0] + cut_span[0]
                    
                cut_site_present = True
                
                return (seg5[0], trim_site, seg5[2]), (trim_site, seg3[1], seg3[2]), cut_site_present

    # There is still overlap, but a cut site could not be detected
    if overlap > 0:
        #print(seg5, seg3)
        if r5_cover:
            trim_site = seg5[1]
        else:
            trim_site = seg3[0]
        return (seg5[0], trim_site, seg5[2]), (trim_site, seg3[1], seg3[2]), cut_site_present


    return seg5, seg3, cut_site_present
    

def trim_mate(seg_keys, original_sequence, ordered_reads, mate, header):

    cut_site_keys = {}
    illegal_post_trim = 0

    if len(seg_keys) == 1:
        existing_data = ordered_reads[0]
        existing_data["adjusted_span"] = seg_keys[0]
        existing_data["trimmed_read"] = existing_data["read"]
        existing_data["new_pairs"] = existing_data["pairs"]
        return ordered_reads, cut_site_keys, illegal_post_trim

    # Chimeric read
    adjusted_seg_keys = seg_keys.copy()
    
    for i in range(len(adjusted_seg_keys)-1):
        seg5 = adjusted_seg_keys[i]
        pairs5 = ordered_reads[i]["pairs"]
        read5 = ordered_reads[i]["read"]
        
        seg3 = adjusted_seg_keys[i+1]
        pairs3 = ordered_reads[i+1]["pairs"]
        read3 = ordered_reads[i+1]["read"]

        mate = ordered_reads[1]["mate"]

        aseg5, aseg3, is_cut = adjust_split(pairs5, seg5, read5,
                                            pairs3, seg3, read3,
                                            original_sequence, mate)

        adjusted_seg_keys[i] = aseg5
        adjusted_seg_keys[i+1] = aseg3
        cut_site_keys[(i, i+1)] = is_cut

    all_keys = list(ordered_reads.keys())
    
    for key in all_keys:
        ordered_reads[key]["adjusted_span"] = adjusted_seg_keys[key]
        
        trimmed_read, new_pairs = create_trimmed_mate(
            ordered_reads[key]["read"], 
            ordered_reads[key]["span"], 
            ordered_reads[key]["adjusted_span"], 
            ordered_reads[key]["pairs"], 
            mate, 
            header
        )

        if trimmed_read == None: # Read is eliminated/erroneous by trimming
            del ordered_reads[key] # Delete read from dictionary 
            illegal_post_trim += 1
        else:
            ordered_reads[key]["trimmed_read"] = trimmed_read
            ordered_reads[key]["new_pairs"] = new_pairs

    return ordered_reads, cut_site_keys, illegal_post_trim
    

def get_5_cut_point(pairs, original_point, is_forward):
    # Handles issues where original cut point is not mapped
    if is_forward:
        index_5 = 0
        iterator = range(len(pairs))
    else:
        index_5 = 1 # Index will be at end of range
        iterator = reversed(range(len(pairs)))
        
    for i in iterator:
        pair = pairs[i]
        if pair["cigar"] == 0:
            # Adjusted pos should be equal to the original point or closer to the 3' end
            if pair["read_pos"] >= original_point: 
                pos_5 = pair["ref_pos"]
                read_5 = pair["read_pos"]
                index_5 += i
                break
    return pos_5, index_5, read_5

def get_3_cut_point(pairs, original_point, is_forward):
    # Handles issues where original cut point is not mapped
    if is_forward:
        index_3 = 1 # Index will be at end of range
        iterator = reversed(range(len(pairs)))
    else:
        index_3 = 0
        iterator = range(len(pairs))
        
    for i in iterator:
        pair = pairs[i]
        if pair["cigar"] == 0:
            # Adjusted pos should be equal to the original point or closer to the 5' end
            if pair["read_pos"] <= original_point:
                pos_3 = pair["ref_pos"]
                read_3 = pair["read_pos"]
                index_3 += i
                break
    return pos_3, index_3, read_3

def build_cigar(read, pairs, start_index, end_index):
    full_cigar = []

    #print(read)
    if not read.is_secondary:
      mask_val = 4
    else:
      mask_val = 5
        
    mask = True
    for i in range(len(pairs)):
        pair = pairs[i]
        cigar_val = pair["cigar"]
        if i < start_index or i >= end_index:
            if cigar_val in [0, 1]:
                cigar_val = mask_val
            elif cigar_val in [2]:
                cigar_val = "*"
        full_cigar.append(cigar_val)

    if read.is_secondary:
        expanded_left = []
        left_original = read.cigartuples[0]
        if left_original[0] == 5:
            expanded_left = [5] * left_original[1]

        expanded_right = []
        right_original = read.cigartuples[-1]
        if right_original[0] == 5:
            expanded_right = [5] * right_original[1]

        full_cigar = expanded_left + full_cigar + expanded_right

    
    updated_cigartuples = []

    #print(full_cigar)
    full_cigar = [i for i in full_cigar if i != "*"]

    counter = 1
    for i in range(1, len(full_cigar)):
        val = full_cigar[i]
        previous_val = full_cigar[i-1]
        if val != previous_val:
            updated_cigartuples.append((previous_val, counter))
            counter = 0
        counter += 1
    updated_cigartuples.append((val, counter))

    return updated_cigartuples


def dedup_contacts(dup_contacts_path, chimeras_path, dedup_contacts_path, save_raw = True):
    columns = ["read", 
               "chrom1", "pos1", "chrom2", 
               "pos2", "strand1", "strand2", 
               "contact_type", "pair_index", "pair_type"]

    contact_stats = {
        "total_contacts" : 0,
        "contacts_dup_rate" : 0,
        "intra1kb_contacts" : 0,
        "intra10kb_contacts" : 0,
        "inter_contacts" : 0,
        "read_pairs_with_contacts" : 0
    }

    contact_types = {"R1" : 0,
                 "R2" : 0,
                 "R1-2" : 0,
                 "R1-2" : 0,
                 "R1&2" : 0,
                 "RU_mask" : 0,
                 "UU_all" : 0,
                 "UU_mask" : 0,
                 "UR_mask" : 0
                 }

    try:
        df = pd.read_table(dup_contacts_path, header=None)
        df.columns = columns
    except ValueError:
        df = pd.DataFrame(columns=columns)
        
    df_dedup = df.drop_duplicates(["chrom1", "pos1", 
                                   "chrom2", "pos2", 
                                   "strand1", "strand2"])

    df_intra = df_dedup[df_dedup["chrom1"] == df_dedup["chrom2"]]
    df_intra_1kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 1000]
    df_intra_10kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 10000]

    df_inter = df_dedup[df_dedup["chrom1"] != df_dedup["chrom2"]]

    if len(df) > 0:
        contact_stats["contacts_dup_rate"] = dedup_rate = 1 - (len(df_dedup) / len(df))
    else:
        contact_stats["contacts_dup_rate"] = np.nan

    contact_stats["intra1kb_contacts"] = len(df_intra_1kb)
    contact_stats["intra10kb_contacts"] = len(df_intra_10kb)
    contact_stats["inter_contacts"] = len(df_inter)
    
    report_contacts = pd.concat([df_intra_1kb, df_inter])

    contact_stats["read_pairs_with_contacts"] = len(report_contacts["read"].unique())
    
    contact_types_emp = {**dict(report_contacts["contact_type"].value_counts()), 
           **dict(report_contacts["pair_type"].value_counts())}
    for i in contact_types_emp:
        contact_types[i] += contact_types_emp[i]

    contact_stats["total_contacts"] = len(report_contacts)

    report_contacts.to_csv(dedup_contacts_path, header=False, index=False, sep="\t")
    
    if not save_raw:
        subprocess.run(shlex.split(f'rm -f {dup_contacts_path}'), check=True)

    contact_stats.update(contact_types)

    try:
        chimeras = pd.read_table(chimeras_path, header=None)
        chimeras.columns = columns
    except ValueError:
        chimeras = pd.DataFrame(columns=columns)
        
    chimeras_dedup = chimeras.drop_duplicates(["chrom1", "pos1", 
                                   "chrom2", "pos2", 
                                   "strand1", "strand2"])

    contact_stats["non_ligation_contacts"] = len(chimeras_dedup)

    return contact_stats

def make_short(contacts_path, short_path, chrom_sizes):

    chrom_order = {}
    index = 0
    with open(chrom_sizes) as c_in:
        for line in c_in:
            line = line.strip().split()
            chrom = line[0]
            chrom_order[chrom] = index
            index += 1

    with gzip.open(contacts_path, 'rt') as contacts_in, \
        gzip.open(short_path, 'wt') as short_out:

        for line in contacts_in:
            line = line.strip().split("\t")
            
            chrom1 = line[1]
            pos1 = line[2]
            chrom2 = line[3]
            pos2 = line[4]
            
            strand1 = line[5]
            strand1 = "0" if strand1 == "+" else "1"
            
            strand2 = line[6]
            strand2 = "0" if strand2 == "+" else "1"

            if chrom1 == chrom2:
                if int(pos1) > int(pos2):
                    short_contact = [strand2, chrom2, pos2, "0", strand1, chrom1, pos1, "1"]
                else:
                    short_contact = [strand1, chrom1, pos1, "0", strand2, chrom2, pos2, "1"]
            elif chrom_order[chrom1] > chrom_order[chrom2]:
                short_contact = [strand2, chrom2, pos2, "0", strand1, chrom1, pos1, "1"]
            else:
                short_contact = [strand1, chrom1, pos1, "0", strand2, chrom2, pos2, "1"]

            short_contact = "\t".join(short_contact) + "\n"

            short_out.write(short_contact)
            

def compute_genome_coverage(bam, min_mapq, min_base_quality, keep_dup=False):
    # Compute genome coverage

    ff_val = "QCFAIL,UNMAP"
    if not keep_dup:
        ff_val = "DUP," + ff_val

    command1 = f"samtools coverage -q {min_mapq} -Q {min_base_quality} -H --ff {ff_val} {bam}"
    command2 = "awk -v OFS='\t' -v nn=0 '{{nn += $5}} END {print nn } '"
    p1 = subprocess.Popen(shlex.split(command1), stdout=subprocess.PIPE)
    p2 = subprocess.run(shlex.split(command2), stdin=p1.stdout, capture_output=True, text=True)
    out = p2.stdout.strip()
    reference_coverage = int(out)
    return reference_coverage

def compute_mapped_nucleotides(bam, min_mapq, min_base_quality, keep_dup=False):
    # Compute mapped nucleotides

    add_flags = "SECONDARY"
    if keep_dup:
        add_flags += ",DUP"
    
    command1 = f"samtools depth -Q {min_mapq} -q {min_base_quality} -g {add_flags} {bam}"
    command2 = "awk -v nn=0 '{nn+=$3} END { print nn}'"
    p1 = subprocess.Popen(shlex.split(command1), stdout=subprocess.PIPE)
    p2 = subprocess.run(shlex.split(command2), stdin=p1.stdout, capture_output=True, text=True)
    out = p2.stdout.strip()
    mapped_read_bases = int(out)
    return mapped_read_bases

def parse_stats(cell, 
                out_prefix, 
                min_mapq = 30, 
                min_base_quality = 20):
    
    trim_stats = out_prefix + "_trim_stats.txt"
    with open(trim_stats) as f:
        line_count = 0
        for line in f:
            if line_count == 1:
                tstats = line.strip().split("\t")
                in_pairs = float(tstats[1])
            elif line_count == 3:
                tstats = line.strip().split("\t")
                out_r1 = float(tstats[6])
            elif line_count == 5:
                tstats = line.strip().split("\t")
                out_r2 = float(tstats[6])
            line_count += 1
    trim_df = pd.DataFrame.from_dict({"demultiplexed_pairs": in_pairs,
                                         "trimmed_R1_mates" : out_r1,
                                         "trimmed_R2_mates" : out_r2
                                        }, orient="index").T

    contam_stats = f"{out_prefix}_contam_stats.txt"
    contam_df = pd.read_table(contam_stats).astype(float)

    dup_stats = f"{out_prefix}_dupsifter_stats.txt"
    dup_count = 0
    with open(dup_stats) as f:
        line_count = 0
        for line in f:
            if line_count in [6, 7]:
                dup_count += int(line.strip().split(": ")[-1])
            line_count += 1
    pre_dedup_count = (contam_df["R1_contam_pass"] + contam_df["R2_contam_pass"]).values[0]
    if pre_dedup_count == 0:
        dup_rate = np.nan
    else:
        dup_rate = dup_count / pre_dedup_count

    dedup_df = pd.DataFrame.from_dict({
        "pre_dedup_mates" : pre_dedup_count,
        "duplicate_mates" : dup_count,
        "alignment_dup_rate" : dup_rate
    }, orient="index").T
                

    alignment_stats = f"{out_prefix}_contacts_trim_stats.txt"
    alignment_df = pd.read_table(alignment_stats)

    methylation_stats = out_prefix + ".allc.tsv.gz_methylation_stats.txt"
    methylation_df = pd.read_table(methylation_stats)

    cell_df = pd.DataFrame([cell], columns=["cell"])

    trimmed_bam = f"{out_prefix}_trimmed_sorted.bam"
    mkdup_bam = f"{out_prefix}_mkdup_sorted.bam"
    masked_bam = f"{out_prefix}_masked_sorted.bam"

    genome_cov_dup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
    genome_cov_dedup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
    genome_cov_dedup_trim = compute_genome_coverage(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
    genome_cov_dedup_mask = compute_genome_coverage(masked_bam, min_mapq, min_base_quality, keep_dup=False)

    mapped_nuc_dup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
    mapped_nuc_dedup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
    mapped_nuc_dedup_trim = compute_mapped_nucleotides(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
    mapped_nuc_dedup_mask = compute_mapped_nucleotides(masked_bam, min_mapq, min_base_quality, keep_dup=False)
    
    cov_map_df = pd.DataFrame([[genome_cov_dup, 
                                genome_cov_dedup,
                                genome_cov_dedup_trim,
                                genome_cov_dedup_mask,
                                mapped_nuc_dup,
                                mapped_nuc_dedup,
                                mapped_nuc_dedup_trim,
                                mapped_nuc_dedup_mask
                               ]], 
                              columns=["reference_coverage_dup",
                                       "reference_coverage_dedup",
                                       "reference_coverage_dedup_trim",
                                       "reference_coverage_dedup_mask",
                                       "mapped_nucleotides_dup",
                                       "mapped_nucleotides_dedup",
                                       "mapped_nucleotides_dedup_trim",
                                       "mapped_nucleotides_dedup_mask",
                                      ]
                             )

    all_stats = pd.concat([cell_df, trim_df, contam_df, dedup_df, alignment_df, cov_map_df, methylation_df], axis=1)

    all_stats = all_stats.T

    out_file = f"{out_prefix}_qc_stats.txt"

    all_stats.to_csv(out_file, sep="\t", header=False)

    