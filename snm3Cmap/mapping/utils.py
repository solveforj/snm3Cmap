import re
import pysam
import pandas as pd
import numpy as np
import subprocess
import shlex
import gzip
import bisect
import os
from collections import OrderedDict

rng = np.random.default_rng(1)

pairs_columns = ["read", 
                 "chrom1", "pos1", "chrom2", 
                 "pos2", "strand1", "strand2", 
                 "pair_type", "rule", "reads", 
                 "contact_class", "overlap", "cut_site_locs"]

def process_chrom_sizes(chrom_sizes_file):
    chrom_sizes = {}
    with open(chrom_sizes_file) as csf:
        for line in csf:
            line = line.strip().split()
            chrom = line[0]
            size = int(line[1])
            chrom_sizes[chrom] = size
    return chrom_sizes
    
def initialize_pairs_file(handle, chrom_sizes):

    handle.write("## pairs format v1.0\n")
    for chrom in chrom_sizes:
        handle.write(f'#chromsize: {chrom} {str(chrom_sizes[chrom])}\n')
    handle.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type rule reads contact_class multimap_overlap cut_site_locs\n")

def make_short(out_prefix):
    
    short_path = f"{out_prefix}_contacts.short.gz"

    try:
        contacts = pd.read_table(pairtools_contacts, comment="#", header=None)
        contacts.columns = pairs_columns
    except:
        contacts = pd.DataFrame(columns = pairs_columns)

    contacts["frag1"] = ["0"] * len(contacts)
    contacts["frag2"] = ["1"] * len(contacts)
    contacts["strand1"] = contacts["strand1"].apply(lambda x : "0" if x == "+" else "1")
    contacts["strand2"] = contacts["strand2"].apply(lambda x : "0" if x == "+" else "1")

    short = contacts[["strand1", "chrom1", "pos1", "frag1", "strand2", "chrom2", "pos2", "frag2"]]
    
    short.to_csv(short_path, header=False, index=False, sep="\t")

def process_restriction_sites(restriction_sites):

    restriction_sites_dict = {}

    for file in restriction_sites:
        enzyme = os.path.basename(file).split("_")[-1].split(".txt")[0]
        restriction_sites_dict[enzyme] = {}
        with open(file) as f:
            for line in f:
                line = line.strip().split()
                chrom = line[0]
                sites = [int(i) for i in line[1:]]
                if chrom not in restriction_sites_dict[enzyme]:
                    restriction_sites_dict[enzyme][chrom] = sites
                else:
                    restriction_sites_dict[enzyme][chrom] += sites

    return restriction_sites_dict

def compute_pair_stats(pair_file, pair_dup_stats=None):

    contact_stats = {
        "total" : 0,
        "dup_rate" : 0,
        "intra1kb" : 0,
        "intra10kb" : 0,
        "intra20kb" : 0,
        "inter" : 0,
        "represented_read_pairs" : 0
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
        contacts = pd.read_table(pair_file, comment="#", header=None)
        contacts.columns = pairs_columns
    except ValueError:
        contacts = pd.DataFrame(columns = pairs_columns)
                    

    contacts["type_rule"] = contacts["pair_type"] + "_" + contacts["rule"]

    df_intra = contacts[contacts["chrom1"] == contacts["chrom2"]]
    df_intra_1kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 1000]
    df_intra_10kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 10000]
    df_intra_20kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 20000]
    df_inter = contacts[contacts["chrom1"] != contacts["chrom2"]]

    contact_stats["intra1kb"] = len(df_intra_1kb)
    contact_stats["intra10kb"] = len(df_intra_10kb)
    contact_stats["intra20kb"] = len(df_intra_20kb)
    contact_stats["inter"] = len(df_inter)
    contact_stats["represented_read_pairs"] = len(contacts["read"].unique())
    
    contact_types_emp = {**dict(contacts["type_rule"].value_counts()), 
           **dict(contacts["reads"].value_counts())}
    for i in contact_types_emp:
        contact_types[i] += contact_types_emp[i]

    for i in contact_types:
        contact_stats[i] = contact_types[i]
    
    contact_stats["total"] = len(contacts)

    pairs_dup_rate = np.nan
    if pair_dup_stats:
        if os.path.isfile(pair_dup_stats):
            with open(pair_dup_stats) as f:
                for line in f:
                    if "summary/frac_dups" in line:
                        try:
                            pairs_dup_rate = float(line.strip().split()[1])
                        except ValueError:
                            pairs_dup_rate = np.nan
                        break
                        
            contact_stats["dup_rate"] = pairs_dup_rate
    
    return contact_stats

def pairtools_stats(out_prefix,
                    contacts,
                    artefacts,
                    contacts_stats=None,
                    artefacts_stats=None,
                    filterbycov_stats=None
                   ):

    stats_path = f"{out_prefix}_pairtools_stats.txt"
    
    contacts_stats = compute_pair_stats(contacts, contacts_stats)

    artefacts_stats = compute_pair_stats(artefacts, artefacts_stats)

    if filterbycov_stats:
        if os.path.isfile(filterbycov_stats):
            with open(filterbycov_stats) as f:
                highcov = 0.0
                for line in f:
                    if "pair_types/FF" in line:
                        try:
                            highcov = float(line.strip().split()[1])
                        except ValueError:
                            highcov = 0.0
                        break
            contacts_stats["high_coverage_pairs"] = highcov
                
    full_stats = {}

    for i in contacts_stats:
        full_stats["contacts_" + i] = contacts_stats[i]

    for i in artefacts_stats:
        full_stats["artefacts_" + i] = artefacts_stats[i]

    stats_df = pd.DataFrame.from_dict(full_stats, orient="index").T
    stats_df.to_csv(stats_path, index=False, sep="\t")
                    
def illegal_overlap(seg_keys):
    illegal_overlap = False
    # If split, determine if sequential overlaps condition is violated
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

def remove_illegal_overlap(seg_keys):
    seg_keys_copy = seg_keys.copy()
    rng.shuffle(seg_keys_copy)
    sorted_seg_keys = sorted(seg_keys_copy, key = lambda x : (x[2], x[1] - x[0]), reverse=True) 

    # This enables skipping already deleted segments
    to_analyze = {}
    for i in sorted_seg_keys:
        to_analyze[i] = True

    for i in range(len(sorted_seg_keys) - 1):
        seg0 = sorted_seg_keys[i]
        if not to_analyze[seg0]:
            continue
            
        for j in range(i + 1, len(sorted_seg_keys)):
            seg1 = sorted_seg_keys[j]

            if seg0[0] >= seg1[0] and seg0[1] <= seg1[1]:
                to_analyze[seg1] = False
            elif seg0[0] <= seg1[0] and seg0[1] >= seg1[1]:
                to_analyze[seg1] = False

    filtered_seg_keys = []
    removed_keys = 0
    for i in seg_keys:
        if to_analyze[i]:
            filtered_seg_keys.append(i)
        else:
            removed_keys += 1

    return filtered_seg_keys, removed_keys

def classify_contact(algn1, 
                     algn2, 
                     pair_index,
                     R1_trimmer,
                     R2_trimmer,
                     rule,
                     restriction_sites,
                     min_inward_dist=1000,
                     min_outward_dist=1000,
                     min_same_strand_dist=0,
                     max_cut_site_whole_algn_dist = 500
                  ):

    R1_cs_keys = R1_trimmer.cut_site_keys
    R2_cs_keys = R2_trimmer.cut_site_keys

    R1_cs_classes = R1_trimmer.cut_site_classes
    R2_cs_classes = R2_trimmer.cut_site_classes

    R1_overlap_keys = R1_trimmer.overlap_keys
    R2_overlap_keys = R2_trimmer.overlap_keys
    
    overlap = 0
    cs_locs = ["na"]
    

    if not algn1["is_mapped"] or not algn1["is_unique"]:
        ct = "na"
        return ct, overlap, cs_locs
    if not algn2["is_mapped"] or not algn2["is_unique"]:
        ct = "na"
        return ct, overlap, cs_locs

    contact_type = pair_index[1]
    contact_class = pair_index[0]

    if rule == "all":
        # Pairtools reports 5' fragment before 3' fragment
        if contact_type == "R1": 
            idx5 = algn1["idx"]
            idx3 = algn2["idx"]
            cs_key = (idx5, idx3)
            if cs_key not in R1_cs_classes:
                ct = "artefact_chimera"
            elif R1_cs_classes[cs_key] == "artefact":
                ct = "artefact_chimera"
                overlap = R1_overlap_keys[cs_key]
                cs_locs = R1_cs_keys[cs_key]
            else:
                ct = R1_cs_classes[cs_key]
                overlap = R1_overlap_keys[cs_key]
                cs_locs = R1_cs_keys[cs_key]
        # Pairtools reports 3' fragment before 5' fragment
        elif contact_type == "R2": 
            idx5 = algn2["idx"]
            idx3 = algn1["idx"]
            cs_key = (idx5, idx3)
            
            if cs_key not in R2_cs_classes:
                ct = "artefact_chimera"
            elif R2_cs_classes[cs_key] == "artefact":
                ct = "artefact_chimera"
                overlap = R2_overlap_keys[cs_key]
                cs_locs = R2_cs_keys[cs_key]
            else:
                ct = R2_cs_classes[cs_key]
                overlap = R2_overlap_keys[cs_key]
                cs_locs = R2_cs_keys[cs_key]
        elif contact_type == "R1&2" or contact_type == "R1-2":
            bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], algn2["read"], 
                                                             restriction_sites, max_cut_site_whole_algn_dist)
            if bp_enzyme != "artefact":
                ct = "gap"
            else:
                ct = "artefact_gap"
                
    elif algn1["type"] == "R":
        if (0, 1) not in R2_cs_classes:
            ct = "artefact_chimera"
        elif R2_cs_classes[(0, 1)] == "artefact":
            ct = "artefact_chimera"
            overlap = R2_overlap_keys[(0, 1)]
            cs_locs = R2_cs_keys[(0, 1)]
        else:
            ct = R2_cs_classes[(0, 1)]
            overlap = R2_overlap_keys[(0, 1)]
            cs_locs = R2_cs_keys[(0, 1)]
    elif algn2["type"] == "R":
        if (0, 1) not in R1_cs_classes:
            ct = "artefact_chimera"
        elif R1_cs_classes[(0, 1)] == "artefact":
            ct = "artefact_chimera"
            overlap = R1_overlap_keys[(0, 1)]
            cs_locs = R1_cs_keys[(0, 1)]
        else:
            ct = R1_cs_classes[(0, 1)]
            overlap = R1_overlap_keys[(0, 1)]
            cs_locs = R1_cs_keys[(0, 1)]

    elif algn1["type"] == "U" and algn2["type"] == "U":
        bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], algn2["read"], 
                                                                     restriction_sites, max_cut_site_whole_algn_dist)
        if bp_enzyme != "artefact":
            ct = "gap"
        else:
            ct = "artefact_gap"

    if algn1["chrom"] == algn2["chrom"]:
        if algn1["pos"] < algn2["pos"]:
            min_algn = algn1
            max_algn = algn2
        else:
            min_algn = algn2
            max_algn = algn1

        if max_algn["strand"] == min_algn["strand"]:
            if (max_algn["pos"] - min_algn["pos"]) < min_same_strand_dist:
                ct = "short"
        elif min_algn["strand"] == "+":
            if (max_algn["pos"] - min_algn["pos"]) < min_inward_dist:
                ct = "short"
        elif min_algn["strand"] == "-":
            if (max_algn["pos"] - min_algn["pos"]) < min_outward_dist:
                ct = "short"

    return ct, overlap, cs_locs

def divide_reads_default(read_group):
    
    r1 = []
    r1_primary = None
    
    r2 = []
    r2_primary = None

    for read in read_group:
        if read.is_read1:
            if not read.is_secondary:
                r1_primary = read
            r1.append(read)
        elif read.is_read2:
            if not read.is_secondary:
                r2_primary = read
            r2.append(read)
        else:
            if not read.is_secondary:
                r1_primary = read
            r1.append(read)

    return r1, r1_primary, r2, r2_primary


def divide_reads_manual_annotation(read_group):
    
    r1 = []
    r1_primary = None
    
    r2 = []
    r2_primary = None

    for read in read_group:
        if read.query_name.split("_")[1] == "1":
            if not read.is_secondary:
                r1_primary = read
            r1.append(read)
        elif read.query_name.split("_")[1] == "2":
            if not read.is_secondary:
                r2_primary = read
            r2.append(read)

    return r1, r1_primary, r2, r2_primary

def get_biscuit_tags(pairs, mate, is_forward):

    ZC = 0
    ZR = 0
    NM = 0
    possible_conv = 0
    
    for pair in pairs:
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
    
    NM = NM - possible_conv

    return ZC, ZR, NM

def get_bwa_tags(pairs):
    NM = 0
    
    for pair in pairs:
        if pair["cigar"] == 0:
            if pair["ref_nuc"].islower():
                NM += 1
        elif pair["cigar"] in [1, 2]:
            NM += 1

    return NM

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
            }
            new_pairs.append(new_algn)
        elif pairs[i]["cigar"] in [2]:
            continue
        else:
            new_pairs.append(pairs[i])
        
    return new_pairs            

def get_loc(read, original_sequence):
    cigar = read.cigartuples
    if not cigar:
        return (0, 0, read.mapping_quality)

    if read.is_forward:
        clip_5 = cigar[0][1] if cigar[0][0] in [4, 5] else 0
        clip_3 = cigar[-1][1] if cigar[-1][0] in [4, 5] else 0
    else:
        clip_5 = cigar[-1][1] if cigar[-1][0] in [4, 5] else 0
        clip_3 = cigar[0][1] if cigar[0][0] in [4, 5] else 0

    start = clip_5
    end = len(original_sequence) - clip_3

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

def get_5_adj_point(pairs, original_point, is_forward):
    # For a given aligned segment of a read, adjust 5' end of segment in case given end is not mapped
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

def get_3_adj_point(pairs, original_point, is_forward):
    # For a given aligned segment of a read, adjust 3' end of segment in case given end is not mapped
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

def get_ref_pos(pairs, new_pos, is_forward):
    # Determines position of pairs where specific reference position is located
    
    # Handles issues where new position is not mapped
    if is_forward:
        index_5 = 0
        iterator = range(len(pairs))
        for i in iterator:
            pair = pairs[i]
            if pair["cigar"] == 0:
                # Adjusted pos should be equal to the new position or closer to the 3' end
                if pair["ref_pos"] >= new_pos: 
                    pos_5 = pair["ref_pos"]
                    read_5 = pair["read_pos"]
                    index_5 += i
                    break
    else:
        index_5 = 1 # Index will be at end of range
        iterator = reversed(range(len(pairs)))
        for i in iterator:
            pair = pairs[i]
            if pair["cigar"] == 0:
                # Adjusted pos should be equal to the new position or closer to the 3' end
                if pair["ref_pos"] <= new_pos: 
                    pos_5 = pair["ref_pos"]
                    read_5 = pair["read_pos"]
                    index_5 += i
                    break

    return pos_5, read_5

def adjust_read5_cut(pairs5, read5, r5_closest, seg5):
    # For the 5' alignment in sequential alignments on the same read, trim alignment range relative to  
    # original read to remove cut site
    if read5.is_forward:
        if r5_closest < read5.reference_end:
            trim_stop = r5_closest - 1
            _, adj_trim_stop = get_ref_pos(pairs5, trim_stop, read5.is_forward)
            seg5 = (seg5[0], adj_trim_stop + 1, seg5[2])
    else:
        if r5_closest > read5.reference_start - 4:
            trim_stop = r5_closest + 4
            _, adj_trim_stop = get_ref_pos(pairs5, trim_stop, read5.is_forward)
            seg5 = (seg5[0], adj_trim_stop + 1, seg5[2])

    return seg5

def adjust_read3_cut(pairs3, read3, r3_closest, seg3):
    # For the 3' alignment in sequential alignments on the same read, trim alignment range relative to 
    # original read to remove cut site
    if read3.is_forward:
        if r3_closest >= read3.reference_start:
            trim_stop = r3_closest + 4
            _, adj_trim_stop = get_ref_pos(pairs3, trim_stop, read3.is_forward)
            seg3 = (adj_trim_stop, seg3[1], seg3[2])
    else:
        if r3_closest <= read3.reference_end - 4:
            trim_stop = r3_closest - 1
            _, adj_trim_stop = get_ref_pos(pairs3, trim_stop, read3.is_forward)
            seg3 = (adj_trim_stop, seg3[1], seg3[2])
            
    return seg3

def get_cut_site_spans(site, seq):
    spans = [i.span() for i in re.finditer(site, seq)]
    return spans


def closest_restriction_site(chrom, pos, restriction_sites_dict, rule="closest"):

    results = {}
    
    for enzyme in restriction_sites_dict:
        results[enzyme] = {}
        
        restriction_sites_chrom = restriction_sites_dict[enzyme][chrom]

        rs_chrom_len = len(restriction_sites_chrom)
        
        insert = bisect.bisect_left(restriction_sites_chrom, pos)

        if insert <= 0:
            fragment = f"{chrom}_{insert}_{insert}"
            
            upstream_pos = restriction_sites_chrom[insert] - 1
            downstream_pos = restriction_sites_chrom[insert] - 1

        elif insert >= rs_chrom_len:
            fragment = f"{chrom}_{insert-1}_{insert-1}"
            
            upstream_pos = restriction_sites_chrom[insert-1] - 1
            downstream_pos = restriction_sites_chrom[insert-1] - 1

        else:
            fragment = f"{chrom}_{insert-1}_{insert}"
    
            upstream_pos = restriction_sites_chrom[insert-1] - 1
            downstream_pos = restriction_sites_chrom[insert] - 1
    
        upstream_dist = np.abs(upstream_pos - pos)
        downstream_dist = np.abs(downstream_pos - pos)

        if rule == "closest":
            if upstream_dist < downstream_dist:
                dist = upstream_dist
                site_pos = upstream_pos
            else:
                dist = downstream_dist
                site_pos = downstream_pos
        elif rule == "upstream":
            dist = upstream_dist
            site_pos = upstream_pos
        elif rule == "downstream":
            dist = downstream_dist
            site_pos = downstream_pos

        results[enzyme] = {"site_pos":site_pos, "dist":dist, "fragment": fragment}

    return results

def classify_breakpoint(r5_site_info, r3_site_info, max_cut_site_distance):

    contact_classes = []
    total_trim = {}
    contact_enzyme = None

    for enzyme in r5_site_info:
        r5_dist = r5_site_info[enzyme]["dist"]
        r3_dist = r3_site_info[enzyme]["dist"]
        r5_less = r5_dist <= max_cut_site_distance
        r3_less = r3_dist <= max_cut_site_distance

        if r5_less and r3_less:
            contact_classes.append(f"{enzyme}_both")
            total_trim[enzyme] = r5_site_info[enzyme]["dist"] + r5_site_info[enzyme]["dist"] 
        elif r5_less:
            contact_classes.append(f"{enzyme}_5")
        elif r3_less:
            contact_classes.append(f"{enzyme}_3")

    # Assign enzyme to valid contacts
    # If there are multiple enzymes assigned to a valid contact, one will be chosen which will require the least trimming
    if len(total_trim) == 0:
        contact_enzyme = "artefact"
    elif len(total_trim) == 1:
        contact_enzyme = list(total_trim.keys())[0]
    elif len(total_trim) > 1:
        min_trim = min(total_trim.values())
        possible_enzymes = [k for k, v in total_trim.items() if v==min_trim]
        if len(possible_enzymes) == 1:
            #contact_enzyme = list(possible_enzymes.keys())[0]
            contact_enzyme = possible_enzymes[0]
        elif len(possible_enzymes) > 1:
            #contact_enzyme = rng.choice(list(possible_enzymes.keys()))
            contact_enzyme = rng.choice(possible_enzymes)
    return contact_classes, contact_enzyme

def split_pair_to_restriction_site(read5, read3, restriction_sites, max_cut_site_distance):

    if read5.is_forward:
        r5_rs = closest_restriction_site(read5.reference_name, 
                                read5.reference_end,
                                restriction_sites)
    else:
        r5_rs = closest_restriction_site(read5.reference_name, 
                        read5.reference_start,
                        restriction_sites) 

    if read3.is_forward:
        r3_rs = closest_restriction_site(read3.reference_name, 
                                        read3.reference_start,
                                        restriction_sites)
    else:
        r3_rs = closest_restriction_site(read3.reference_name, 
                                        read3.reference_end,
                                        restriction_sites)

    bp_class, bp_enzyme = classify_breakpoint(r5_rs, r3_rs, max_cut_site_distance)

    return bp_class, bp_enzyme, r5_rs, r3_rs

def gap_pair_to_restriction_site(read5, read3, restriction_sites, max_cut_site_distance):

    if read5.is_forward:
        r5_rule = "downstream"
        r5_pos = read5.reference_end
    else:
        r5_rule = "upstream"
        r5_pos = read5.reference_start

    r5_rs = closest_restriction_site(read5.reference_name, 
                                     r5_pos, restriction_sites, rule=r5_rule) 

    if read3.is_forward:
        r3_rule = "upstream"
        r3_pos = read3.reference_start
    else:
        r3_rule = "downstream"
        r3_pos = read3.reference_end
        
    r3_rs = closest_restriction_site(read3.reference_name, 
                                     r3_pos, restriction_sites, rule=r3_rule) 
    
    bp_class, bp_enzyme = classify_breakpoint(r5_rs, r3_rs, max_cut_site_distance)

    return bp_class, bp_enzyme, r5_rs, r3_rs

class ReadTrimmer:

    def adjust_split(self,
                     pairs5, 
                     seg5, 
                     read5,
                     pairs3,
                     seg3,
                     read3
                    ):

        original_sequence = self.original_sequence
        mate = self.mate
        
        overlap = seg5[1] - seg3[0]
        original_overlap = overlap
        cut_site_present = "none"
        
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


        read5_cut = False
        read3_cut = False
        

        bp_class, bp_enzyme, r5_rs, r3_rs = split_pair_to_restriction_site(read5, read3, 
                                                       self.restriction_sites, 
                                                       self.max_cut_site_split_algn_dist)
        self.total_pairs += 1

        if overlap > 0:
            self.first_trim += 1
        
        if bp_enzyme != "artefact":
            r5_closest = r5_rs[bp_enzyme]["site_pos"]
            r3_closest = r3_rs[bp_enzyme]["site_pos"]
            seg5 = adjust_read5_cut(pairs5, read5, r5_closest, seg5)
            seg3 = adjust_read3_cut(pairs3, read3, r3_closest, seg3)
            overlap = seg5[1] - seg3[0]
            if overlap > 0:
                self.second_trim += 1
            self.cut_site_pairs += 1
        
        if overlap > 0:
            # Handle overlap even if both reads had cut site trimming or neither reads had cut site
            if r5_cover:
                # Trim 3' read
                seg3 = (seg5[1], seg3[1], seg3[2])
            else:
                seg5 = (seg5[0], seg3[0], seg5[2])

        return seg5, seg3, bp_class, bp_enzyme, original_overlap
    
    def create_trimmed_mate(self, read, original_span, adjusted_span, pairs):

        mate = self.mate
        header = self.header
        full_bam = self.full_bam
        
        if original_span == adjusted_span:
            return read, pairs
            
        pos_5, index_5, read_5 = get_5_adj_point(pairs, adjusted_span[0], read.is_forward)
        # Index is at end of range; subtract 1 to get actual final index
        pos_3, index_3, read_3 = get_3_adj_point(pairs, adjusted_span[1] - 1, read.is_forward)
    
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
        new_MD = get_md_tag(new_pairs)

        if self.bisulfite:
            new_ZC, new_ZR, new_NM = get_biscuit_tags(new_pairs, mate, read.is_forward)
    
            new_tags = {
                "ZC" : new_ZC,
                "ZR" : new_ZR,
                "NM" : new_NM,
                "MD" : new_MD
            }
        else:
            new_NM = get_bwa_tags(new_pairs)

            new_tags = {
                "NM" : new_NM
            }
    
        tags = read.get_tags()
    
        for i in range(len(tags)):
            tag = tags[i]
            if tag[0] in new_tags:
                tags[i] = (tag[0], new_tags[tag[0]])
    
        real_trim_5 = read_5 - original_span[0]
        real_trim_3 = original_span[1] - read_3 - 1 
    
        if full_bam:
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
    
    def trim_mate(self):

        ordered_reads = self.ordered_reads
        seg_keys = self.seg_keys
        mate = self.mate
        
        cut_site_keys = {}
        overlap_keys = {}
        cut_site_labels = {}
        cut_site_classes = {}
        illegal_post_trim = 0
    
        if len(seg_keys) == 1:
            existing_data = ordered_reads[0]
            existing_data["adjusted_span"] = seg_keys[0]
            existing_data["trimmed_read"] = existing_data["read"]
            existing_data["new_pairs"] = existing_data["pairs"]

            self.cut_site_keys = cut_site_keys
            self.cut_site_labels = cut_site_labels
            self.cut_site_classes = cut_site_classes
            self.overlap_keys = overlap_keys
            self.illegal_post_trim = illegal_post_trim
            
            return
    
        # split read
        adjusted_seg_keys = seg_keys.copy()
        
        for i in range(len(adjusted_seg_keys)-1):
            seg5 = adjusted_seg_keys[i]
            pairs5 = ordered_reads[i]["pairs"]
            read5 = ordered_reads[i]["read"]
            
            seg3 = adjusted_seg_keys[i+1]
            pairs3 = ordered_reads[i+1]["pairs"]
            read3 = ordered_reads[i+1]["read"]
    
            mate = ordered_reads[1]["mate"]
    
            aseg5, aseg3, bp_class, bp_enzyme, overlap = self.adjust_split(pairs5, 
                                                                            seg5, 
                                                                            read5, 
                                                                            pairs3, 
                                                                            seg3, 
                                                                            read3)
    
            adjusted_seg_keys[i] = aseg5
            adjusted_seg_keys[i+1] = aseg3
            cut_site_keys[(i, i+1)] = bp_class
            cut_site_classes[(i, i+1)] = bp_enzyme
            overlap_keys[(i, i+1)] = overlap
    
            if i not in cut_site_labels:
                cut_site_labels[i] = []
            if i+1 not in cut_site_labels:
                cut_site_labels[i+1] = []
    
            if bp_enzyme != "artefact":
                cut_site_labels[i].append("D")
                cut_site_labels[i+1].append("U")
    
        all_keys = list(ordered_reads.keys())
        
        for key in all_keys:
            ordered_reads[key]["adjusted_span"] = adjusted_seg_keys[key]
            
            trimmed_read, new_pairs = self.create_trimmed_mate(
                ordered_reads[key]["read"], 
                ordered_reads[key]["span"], 
                ordered_reads[key]["adjusted_span"], 
                ordered_reads[key]["pairs"]
            )
    
            if trimmed_read == None: # Read is eliminated/erroneous by trimming
                del ordered_reads[key] # Delete read from dictionary 
                illegal_post_trim += 1
            else:
                ordered_reads[key]["trimmed_read"] = trimmed_read
                ordered_reads[key]["new_pairs"] = new_pairs

        self.cut_site_keys = cut_site_keys
        self.cut_site_labels = cut_site_labels
        self.cut_site_classes = cut_site_classes
        self.overlap_keys = overlap_keys
        self.illegal_post_trim = illegal_post_trim
        
        return
        
    def process_mate(self):

        seg_keys = self.seg_keys
        has_split = self.has_split
        read_parts = self.read_parts
        original_sequence = self.original_sequence
        mate = self.mate
        
        ordered_reads = OrderedDict()

        if len(seg_keys) == 1 and not has_split:
            read = read_parts[seg_keys[0]]
            span = seg_keys[0]
            ordered_reads[0] = {"read" : read,
                                "span" : span,
                                "pairs" : None,
                                "mate" : mate,
                                "adjusted_span" : span,
                                "trimmed_read" : read,
                                "new_pairs" : None}

            self.ordered_reads = ordered_reads
            self.cut_site_keys = {}
            self.cut_site_labels = {}
            self.cut_site_classes = {}
            self.overlap_keys = {}
            self.illegal_post_trim = 0

            return
            
        for i in range(len(seg_keys)):
            read = read_parts[seg_keys[i]]
            span = seg_keys[i]
            ordered_reads[i] = {"read" : read,
                                "span" : span,
                                "pairs" : alignment_info(read, span, original_sequence),
                                "mate" : mate
                               }

        self.ordered_reads = ordered_reads

        self.trim_mate()

        return
    
    def __init__(self,
                 seg_keys, 
                 original_sequence = "",
                 read_parts = {},
                 restriction_sites = {},
                 header = None,
                 full_bam = False,
                 mate = "1",
                 bisulfite=False,
                 max_cut_site_split_algn_dist=20,
                 max_cut_site_whole_algn_dist=500
                ):

        if len(seg_keys) == 0:
            self.ordered_reads = OrderedDict()
            self.has_split = False
            self.cut_site_keys = {}
            self.cut_site_labels = {}
            self.overlap_keys = {}
            self.cut_site_classes = {}
            self.original_sequence = None
            self.cut_site_pairs = 0
            self.first_trim = 0
            self.second_trim = 0
            self.total_pairs = 0
            return
            
        self.has_split = read_parts[seg_keys[0]].has_tag("SA")
            
        self.seg_keys = seg_keys
        self.original_sequence = original_sequence
        self.read_parts = read_parts
        self.restriction_sites = restriction_sites
        self.header = header
        self.full_bam = full_bam
        self.mate = mate
        self.bisulfite = bisulfite
        self.max_cut_site_split_algn_dist = max_cut_site_split_algn_dist
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist
        self.cut_site_pairs = 0
        self.first_trim = 0
        self.second_trim = 0
        self.total_pairs = 0
        

        self.process_mate()
            
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

def snm3Cseq_qc_stats(job,
                       out_prefix, 
                       min_mapq = 30, 
                       min_base_quality = 20):
    
    txt_paths = [f"{out_prefix}_trim_stats.txt",
                 f"{out_prefix}_contam_stats.txt",
                 f"{out_prefix}_dupsifter_stats.txt", 
                 f"{out_prefix}_alignment_stats.txt",
                 f"{out_prefix}_pairtools_stats.txt",
                 f"{out_prefix}.allc.tsv.gz_methylation_stats.txt"
                ]
    
    stat_dfs = [pd.DataFrame([job], columns=["job"])]
    for path in txt_paths:
        if os.path.exists(path):
            stat_dfs.append(pd.read_table(path).astype(float))                 
    
    cov_map_stats = []
    cov_map_columns = []
    trimmed_bam = f"{out_prefix}_trimmed_sorted.bam"
    if os.path.exists(trimmed_bam):
        genome_cov_dedup_trim = compute_genome_coverage(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dedup_trim = compute_mapped_nucleotides(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dedup_trim, mapped_nuc_dedup_trim]
        cov_map_columns += ["reference_coverage_dedup_trim", "mapped_nucleotides_dedup_trim"]

    mkdup_bam = f"{out_prefix}_mkdup_sorted.bam"
    if os.path.exists(mkdup_bam):
        genome_cov_dup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
        genome_cov_dedup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
        mapped_nuc_dedup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dup, genome_cov_dedup, mapped_nuc_dup, mapped_nuc_dedup]
        cov_map_columns += ["reference_coverage_dup", "reference_coverage_dedup", "mapped_nucleotides_dup", "mapped_nucleotides_dedup"]

    masked_bam = f"{out_prefix}_masked_sorted.bam"
    if os.path.exists(masked_bam):
        genome_cov_dedup_mask = compute_genome_coverage(masked_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dedup_mask = compute_mapped_nucleotides(masked_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dedup_mask, mapped_nuc_dedup_mask]
        cov_map_columns += ["reference_coverage_dedup_mask", "mapped_nucleotides_dedup_mask"]

    cov_map_df = pd.DataFrame([cov_map_stats], columns=cov_map_columns)
    stat_dfs.append(cov_map_df)

    return stat_dfs


    
def aggregate_qc_stats(job,
                       out_prefix, 
                       mode,
                       min_mapq = 30, 
                       min_base_quality = 20):

    if mode in ["snm3Cseq", "bsdna"]:

        txt_paths = [f"{out_prefix}_trim_stats.txt",
                 f"{out_prefix}_contam_stats.txt",
                 f"{out_prefix}_dupsifter_stats.txt", 
                 f"{out_prefix}_alignment_stats.txt",
                 f"{out_prefix}_pairtools_stats.txt",
                 f"{out_prefix}.allc.tsv.gz_methylation_stats.txt"
                ]

    elif mode == "dna":
        
        txt_paths = [f"{out_prefix}_trim_stats.txt",
                 f"{out_prefix}_alignment_stats.txt",
                 f"{out_prefix}_pairtools_stats.txt",
                ]
    
    stat_dfs = [pd.DataFrame([job], columns=["job"])]
    for path in txt_paths:
        if os.path.exists(path):
            stat_dfs.append(pd.read_table(path).astype(float))                 
    
    cov_map_stats = []
    cov_map_columns = []
    trimmed_bam = f"{out_prefix}_trimmed_sorted.bam"
    if os.path.exists(trimmed_bam):
        genome_cov_dedup_trim = compute_genome_coverage(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dedup_trim = compute_mapped_nucleotides(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dedup_trim, mapped_nuc_dedup_trim]
        cov_map_columns += ["reference_coverage_dedup_trim", "mapped_nucleotides_dedup_trim"]

    mkdup_bam = f"{out_prefix}_mkdup_sorted.bam"
    if os.path.exists(mkdup_bam):
        genome_cov_dup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
        genome_cov_dedup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
        mapped_nuc_dedup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dup, genome_cov_dedup, mapped_nuc_dup, mapped_nuc_dedup]
        cov_map_columns += ["reference_coverage_dup", "reference_coverage_dedup", "mapped_nucleotides_dup", "mapped_nucleotides_dedup"]

    masked_bam = f"{out_prefix}_masked_sorted.bam"
    if os.path.exists(masked_bam):
        genome_cov_dedup_mask = compute_genome_coverage(masked_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dedup_mask = compute_mapped_nucleotides(masked_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dedup_mask, mapped_nuc_dedup_mask]
        cov_map_columns += ["reference_coverage_dedup_mask", "mapped_nucleotides_dedup_mask"]

    cov_map_df = pd.DataFrame([cov_map_stats], columns=cov_map_columns)
    stat_dfs.append(cov_map_df)

    all_stats = pd.concat(stat_dfs, axis=1)
    all_stats = all_stats.T

    out_file = f"{out_prefix}_qc_stats.txt"

    all_stats.to_csv(out_file, sep="\t", header=False)

    