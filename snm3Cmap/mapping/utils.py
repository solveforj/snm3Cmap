import re
import pysam
import pandas as pd
import numpy as np

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
                   pair_index):
    if not algn1["is_mapped"] or not algn1["is_unique"]:
        return False
    if not algn2["is_mapped"] or not algn2["is_unique"]:
        return False
    contact_type = pair_index[1]
    
    if contact_type == "R1": # Pairtools reports 5' fragment before 3' fragment
        if algn1["is_cut"] != True:
            return False
    elif contact_type == "R2": # Pairtools reports 3' fragment before 5' fragment
        if algn2["is_cut"] != True:
            return False
        
    return True

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

    if len(seg_keys) == 1:
        existing_data = ordered_reads[0]
        existing_data["adjusted_span"] = seg_keys[0]
        existing_data["is_cut"] = False
        existing_data["trimmed_read"] = existing_data["read"]
        existing_data["new_pairs"] = existing_data["pairs"]

    # Chimeric read
    adjusted_seg_keys = seg_keys.copy()
    cut_site_keys = [False] * len(seg_keys)
    
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

        #print(aseg5, aseg3)

        adjusted_seg_keys[i] = aseg5
        adjusted_seg_keys[i+1] = aseg3
        cut_site_keys[i] = is_cut

    for i in ordered_reads:
        existing_data = ordered_reads[i]
        existing_data["adjusted_span"] = adjusted_seg_keys[i]
        existing_data["is_cut"] = cut_site_keys[i]
        #for pair in existing_data["pairs"]:
        #    print(pair)
        trimmed_read, new_pairs = create_trimmed_mate(existing_data["read"], 
                                                      existing_data["span"], 
                                                      existing_data["adjusted_span"], 
                                                      existing_data["pairs"], 
                                                      mate, header)
        existing_data["trimmed_read"] = trimmed_read
        existing_data["new_pairs"] = new_pairs
    

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

def parse_stats(cell, out_prefix):
    
    trim_stats = out_prefix + "_trim_stats.txt"
    with open(trim_stats) as f:
        line_count = 0
        for line in f:
            if line_count == 1:
                tstats = line.strip().split("\t")
                in_pairs = tstats[1]
            elif line_count == 3:
                tstats = line.strip().split("\t")
                out_r1 = tstats[6]
            elif line_count == 5:
                tstats = line.strip().split("\t")
                out_r2 = tstats[6]
            line_count += 1
    trim_df = pd.DataFrame.from_dict({"demultiplexed_pairs": in_pairs,
                                         "trimmed_R1_mates" : out_r1,
                                         "trimmed_R2_mates" : out_r2
                                        }, orient="index").T

    contam_stats = out_prefix + "_contam_stats.txt"
    contam_df = pd.read_table(contam_stats)

    dup_stats = out_prefix + "_dupsifter_stats.txt"
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
        "algn_dup_rate" : dup_rate
    }, orient="index").T
                

    alignment_stats = out_prefix + "_alignment_stats.txt"
    alignment_df = pd.read_table(alignment_stats)
    
    chrl_stats = out_prefix + ".allc.tsv.gz_chrL_stats.txt"
    chrl_df = pd.read_table(chrl_stats)

    all_stats = pd.concat([trim_df, contam_df, dedup_df, alignment_df, chrl_df], axis=1)

    all_stats.index = [cell]

    out_file = out_prefix + "_qc_stats.txt"

    all_stats.to_csv(out_file, sep="\t")

    