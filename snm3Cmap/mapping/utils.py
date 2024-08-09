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

def process_chrom_sizes(chrom_sizes_file):
    chrom_sizes = {}
    with open(chrom_sizes_file) as csf:
        for line in csf:
            line = line.strip().split()
            chrom = line[0]
            size = int(line[1])
            chrom_sizes[chrom] = size
    return chrom_sizes

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

def aggregate_qc_stats(job,
                       out_prefix, 
                       mode,
                       min_mapq = 30, 
                       min_base_quality = 20):

    if mode == "bsdna":
        
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

    