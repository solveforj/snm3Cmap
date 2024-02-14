import pysam
import gzip
import numpy as np
from collections import OrderedDict
import subprocess
import shlex
import pandas as pd

from .utils import *
from .pairtools import *

class ContactGenerator:
    
    def process_mate(self, all_alignments, primary_alignment, mate, read_group_name):

        if len(all_alignments) == 0:
            return OrderedDict(), {}, None, 0, 0, 0
    
        read_parts = {}
        
        primary_alignment = primary_alignment
    
        # Get original whole read sequence from primary alignment
        original_sequence = primary_alignment.get_forward_sequence()
    
        for read in all_alignments:
            if read.is_secondary:
                if "S" in read.cigarstring:
                    read.flag = read.flag - 256
            read_parts[get_loc(read)] = read
    
        # Sort segments in order from 5' to 3' by starting position
        # Only choose segments with high MAPQ
        seg_keys = sorted([i for i in list(read_parts.keys()) if i[2] >= self.min_mapq], 
                          key = lambda x : x[0])


        if len(seg_keys) == 0:
            return OrderedDict(), {}, None, 0, 0, 0

        total_alignments = len(seg_keys)
        has_chimera = read_parts[seg_keys[0]].has_tag("SA")
        whole_alignments = 0 if has_chimera else 1
        chimeric_alignments = 1 if has_chimera else 0
        
        # Throw out these reads

        if illegal_overlap(seg_keys):
            self.illegal_overlap_count += 1
            return OrderedDict(), {}, None, total_alignments, whole_alignments, chimeric_alignments

        ordered_reads = OrderedDict()

        if len(seg_keys) == 1 and not has_chimera:
            read = read_parts[seg_keys[0]]
            span = seg_keys[0]
            ordered_reads[0] = {"read" : read,
                                "span" : span,
                                "pairs" : None,
                                "mate" : mate,
                                "adjusted_span" : span,
                                "trimmed_read" : read,
                                "new_pairs" : None}
            cut_site_keys = {}
        else:
            for i in range(len(seg_keys)):
                read = read_parts[seg_keys[i]]
                span = seg_keys[i]
                ordered_reads[i] = {"read" : read,
                                    "span" : span,
                                    "pairs" : alignment_info(read, span, original_sequence),
                                    "mate" : mate
                                   }
                
            ordered_reads, cut_site_keys = trim_mate(seg_keys, original_sequence, ordered_reads, mate, self.header)
        
        return ordered_reads, cut_site_keys, original_sequence, total_alignments, whole_alignments, chimeric_alignments
    
    def process_read_group(self, read_group, read_group_name, bam_out, contacts_out):
        r1 = []
        r1_primary = None
        
        r2 = []
        r2_primary = None

        # Divide reads into R1 and R2
        for read in read_group:
            if read.query_name.split("_")[1] == "1":
                if not read.is_secondary:
                    r1_primary = read
                r1.append(read)
            elif read.query_name.split("_")[1] == "2":
                if not read.is_secondary:
                    r2_primary = read
                r2.append(read)

        # Order reads from 5' to 3'
        R1, R1_cs_keys, R1_oseq, R1_tot_algn, R1_w_algn, R1_c_algn = self.process_mate(r1, r1_primary, "1", read_group_name)
        R2, R2_cs_keys, R2_oseq, R2_tot_algn, R2_w_algn, R2_c_algn = self.process_mate(r2, r2_primary, "2", read_group_name)

        self.r1_total_alignments += R1_tot_algn
        self.r1_whole_alignments += R1_w_algn
        self.r1_chimeric_alignments += R1_c_algn
        
        self.r2_total_alignments += R2_tot_algn
        self.r2_whole_alignments += R2_w_algn
        self.r2_chimeric_alignments += R2_c_algn
        
        contacts, rule = contact_iter(R1, R2, min_mapq=30, max_molecule_size=750, max_inter_align_gap=20)

        for (hic_algn1, hic_algn2, pair_index) in contacts:
            if contact_filter(hic_algn1, hic_algn2, pair_index, R1_cs_keys, R2_cs_keys):
                write_pairsam(hic_algn1, hic_algn2, read_group_name, pair_index, rule, contacts_out)

        for i in R1:
            read = R1[i]["trimmed_read"]
            if not read.is_duplicate:
                bam_out.write(read)
            else:
                self.r1_duplicate_alignments += 1

        for i in R2:
            read = R2[i]["trimmed_read"]
            if not read.is_duplicate:
                bam_out.write(read)
            else:
                self.r2_duplicate_alignments += 1
   
    def process_bam(self):

        self.illegal_overlap_count = 0

        self.r1_total_alignments = 0
        self.r1_whole_alignments = 0
        self.r1_chimeric_alignments = 0
        
        self.r2_total_alignments = 0
        self.r2_whole_alignments = 0
        self.r2_chimeric_alignments = 0

        self.r1_duplicate_alignments = 0
        self.r2_duplicate_alignments = 0
        
        iter_count = 0
        count = 0
        with pysam.AlignmentFile(self.bam, index_filename=None) as bam_in:
            with pysam.AlignmentFile(self.masked_bam, 'wb', template=bam_in) as bam_out, \
                gzip.open(self.contacts, 'wt') as contacts_out:

                self.header = bam_in.header.to_dict()
                
                for read in bam_in:
                    read_id = read.query_name.split("_")[0]
                    if iter_count == 0:
                        read_group_name = read_id
                        read_group = []
                        iter_count = 1
                    count += 1
                    
                    if read_group_name == read_id:
                        read_group.append(read)
                    else:
                        self.process_read_group(read_group, read_group_name, bam_out, contacts_out)
                        
                        read_group_name = read_id
                        read_group = [read]
                self.process_read_group(read_group, read_group_name, bam_out, contacts_out)

    def dedup_contacts(self, save_raw = True):
        columns = ["read", 
                   "chrom1", "pos1", "chrom2", 
                   "pos2", "strand1", "strand2", 
                   "contact_type", "pair_index", "pair_type"]

        try:
            df = pd.read_table(self.contacts, header=None)
            df.columns = columns
        except ValueError:
            df = pd.DataFrame(columns=columns)
            
        df_dedup = df.drop_duplicates(["chrom1", "pos1", 
                                       "chrom2", "pos2", 
                                       "strand1", "strand2"])

        if len(df) > 0:
            self.dedup_rate = 1 - (len(df_dedup) / len(df))
        else:
            self.dedup_rate = np.nan

        df_intra = df_dedup[df_dedup["chrom1"] == df_dedup["chrom2"]]
        df_intra_1kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 1000]
        df_intra_10kb = df_intra[(df_intra["pos2"] - df_intra["pos1"]).abs() >= 10000]

        df_inter = df_dedup[df_dedup["chrom1"] != df_dedup["chrom2"]]
        
        self.intra_1kb = len(df_intra_1kb)
        self.intra_10kb = len(df_intra_10kb)
        self.inter = len(df_inter)

        self.contact_types = {"R1" : 0,
                              "R2" : 0,
                              "R1-2" : 0,
                              "R1-2" : 0,
                              "R1&2" : 0,
                              "RU_mask" : 0,
                              "UU_all" : 0,
                              "UU_mask" : 0,
                              "UR_mask" : 0
                             }

        report_contacts = pd.concat([df_intra_1kb, df_inter])
        contacts_stats = {**dict(report_contacts["contact_type"].value_counts()), 
               **dict(report_contacts["pair_type"].value_counts())}
        for i in contacts_stats:
            self.contact_types[i] += contacts_stats[i]

        self.total_contacts = len(report_contacts)

        report_contacts.to_csv(self.contacts_dedup, header=False, index=False, sep="\t")
        
        if not save_raw:
            subprocess.run(shlex.split(f'rm -f {self.contacts}'), check=True)

    def generate_stats(self):

        total_alignments = self.r1_total_alignments + self.r2_total_alignments
        duplicated_alignments = self.r1_duplicate_alignments + self.r2_duplicate_alignments

        if total_alignments > 0:
            alignments_dedup_rate = duplicated_alignments / total_alignments 
        else:
            alignments_dedup_rate = np.nan
            
        stats_dict = {
            "R1_total_alignments" : self.r1_total_alignments,
            #"R1_duplicate_alignments" : self.r1_duplicate_alignments,
            "R1_whole_algn_mates" : self.r1_whole_alignments,
            "R1_chimeric_algn_mates" : self.r1_chimeric_alignments,
            "R1_total_algn_mates" : self.r1_chimeric_alignments +  self.r1_whole_alignments,
            
            "R2_total_alignments" : self.r2_total_alignments,
            #"R2_duplicate_alignments" : self.r2_duplicate_alignments,
            "R2_whole_algn_mates" : self.r2_whole_alignments,
            "R2_chimeric_algn_mates" : self.r2_chimeric_alignments,
            "R2_total_algn_mates" : self.r2_chimeric_alignments +  self.r2_whole_alignments,

            "discarded_chimeric_mates" : self.illegal_overlap_count,
            #"alignments_dedup_rate" : alignments_dedup_rate,

            "total_contacts" : self.total_contacts,
            "contacts_dedup_rate" : self.dedup_rate,
            "intra1kb_contacts" : self.intra_1kb,
            "intra10kb_contacts" : self.intra_10kb,
            "inter_contacts" : self.inter
        }

        stats_dict.update(self.contact_types)
        
        stats_df = pd.DataFrame.from_dict(stats_dict, orient="index").T

        stats_df.to_csv(self.stats, index=False, sep="\t")
        
    def __init__(self, bam, out_prefix):

        self.min_mapq = 30
        self.max_molecule_size = 750
        self.max_inter_align_gap = 20

        self.contacts = f'{out_prefix}_contacts.txt.gz'
        self.contacts_dedup = f'{out_prefix}_contacts_dedup.txt.gz'
        self.stats = f"{out_prefix}_alignment_stats.txt"
        
        self.bam = bam
        self.masked_bam = f'{out_prefix}_trimmed.bam'

        self.process_bam()

        self.dedup_contacts(save_raw=False)
        
        self.generate_stats()