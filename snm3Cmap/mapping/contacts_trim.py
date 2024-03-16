import pysam
import gzip
import numpy as np
from collections import OrderedDict

from .utils import *
from .pairtools import *

class ContactGenerator:
    
    def process_mate(self, all_alignments, primary_alignment, mate, read_group_name):

        if len(all_alignments) == 0:
            return OrderedDict(), {}, None, False#, 0, 0, 0
    
        read_parts = {}
            
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
            return OrderedDict(), {}, None, False#, 0, 0, 0

        has_chimera = read_parts[seg_keys[0]].has_tag("SA")

        # Throw out these reads
        if illegal_overlap(seg_keys):
            self.illegal_overlap_count += len(seg_keys)
            return OrderedDict(), {}, None, has_chimera#, total_alignments, whole_alignments, chimeric_alignments

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
                
            ordered_reads, cut_site_keys, illegal_post_trim = trim_mate(seg_keys, original_sequence, ordered_reads, mate, self.header)
            self.illegal_post_trim_count += illegal_post_trim
        
        return ordered_reads, cut_site_keys, original_sequence, has_chimera#, total_alignments, whole_alignments, chimeric_alignments
    
    def process_read_group(self, read_group, read_group_name, bam_out, contacts_out, chimeras_out):
        r1 = []
        r1_primary = None
        
        r2 = []
        r2_primary = None

        # Divide reads into R1 and R2
        for read in read_group:
            if read.query_name.split("_")[1] == "1":
                if not read.is_secondary:
                    r1_primary = read
                    r1_primary_unique = read.mapping_quality >= self.min_mapq
                r1.append(read)
            elif read.query_name.split("_")[1] == "2":
                if not read.is_secondary:
                    r2_primary = read
                    r2_primary_unique = read.mapping_quality >= self.min_mapq
                r2.append(read)

        # Order reads from 5' to 3'
        R1, R1_cs_keys, R1_oseq, R1_has_chimera = self.process_mate(r1, r1_primary, "1", read_group_name)
        R2, R2_cs_keys, R2_oseq, R2_has_chimera = self.process_mate(r2, r2_primary, "2", read_group_name)

        if len(R1) + len(R2) >= 2:
            self.at_least_two_alignments += 1
        
        contacts, rule = contact_iter(R1, R2, min_mapq=30, max_molecule_size=750, max_inter_align_gap=20)

        for (hic_algn1, hic_algn2, pair_index) in contacts:
            contact_class = contact_filter(hic_algn1, hic_algn2, pair_index, R1_cs_keys, R2_cs_keys)
            if contact_class == "contact":
                write_pairsam(hic_algn1, hic_algn2, read_group_name, pair_index, rule, contacts_out)
            elif contact_class == "chimera":
                write_pairsam(hic_algn1, hic_algn2, read_group_name, pair_index, rule, chimeras_out)

        # dupsifter marks all secondary alignments as duplicate if primary alignment is duplicate
        
        # Write R1 reads
        if len(R1) > 0:
            if not r1_primary.is_duplicate:
                if R1_has_chimera:
                    self.r1_chimeric_alignments_dedup += 1
                else:
                    self.r1_whole_alignments_dedup += 1
    
                # Read is chimeric, but primary alignment will not be included in processed mates
                if R1_has_chimera and not r1_primary_unique:
                    bam_out.write(r1_primary)

            if R1_has_chimera:
                self.r1_chimeric_alignments_dup += 1
            else:
                self.r1_whole_alignments_dup += 1

            for i in R1:
                read = R1[i]["trimmed_read"]
                if not read.is_duplicate:
                    self.r1_total_alignments_dedup += 1
                    bam_out.write(read)
                
                self.r1_total_alignments_dup += 1

        # Write R2 reads
        if len(R2) > 0:
            if not r2_primary.is_duplicate:
                if R2_has_chimera:
                    self.r2_chimeric_alignments_dedup += 1
                else:
                    self.r2_whole_alignments_dedup += 1
    
                # Read is chimeric, but primary alignment will not be included in processed mates
                if R2_has_chimera and not r2_primary_unique:
                    bam_out.write(r2_primary)

            if R2_has_chimera:
                self.r2_chimeric_alignments_dup += 1
            else:
                self.r2_whole_alignments_dup += 1

            for i in R2:
                read = R2[i]["trimmed_read"]
                if not read.is_duplicate:
                    self.r2_total_alignments_dedup += 1
                    bam_out.write(read)

                self.r2_total_alignments_dup += 1
   
    def process_bam(self):

        self.illegal_overlap_count = 0
        self.illegal_post_trim_count = 0
        self.at_least_two_alignments = 0
        
        self.r1_total_alignments_dup = 0
        self.r1_whole_alignments_dup = 0
        self.r1_chimeric_alignments_dup = 0

        self.r1_total_alignments_dedup = 0
        self.r1_whole_alignments_dedup = 0
        self.r1_chimeric_alignments_dedup = 0
        
        self.r2_total_alignments_dup = 0
        self.r2_whole_alignments_dup = 0
        self.r2_chimeric_alignments_dup = 0

        self.r2_total_alignments_dedup = 0
        self.r2_whole_alignments_dedup = 0
        self.r2_chimeric_alignments_dedup = 0
        
        iter_count = 0
        count = 0
        with pysam.AlignmentFile(self.bam, index_filename=None) as bam_in:
            with pysam.AlignmentFile(self.masked_bam, 'wb', template=bam_in) as bam_out, \
                gzip.open(self.contacts, 'wt') as contacts_out, \
                gzip.open(self.chimeras, 'wt') as chimeras_out:

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
                        self.process_read_group(read_group, read_group_name, 
                                                bam_out, contacts_out, chimeras_out)
                        
                        read_group_name = read_id
                        read_group = [read]
                self.process_read_group(read_group, read_group_name, 
                                        bam_out, contacts_out, chimeras_out)

    def generate_stats(self):

        stats_dict = {
            "R1_total_alignments_dup" : self.r1_total_alignments_dup,
            "R1_whole_aligned_mates_dup" : self.r1_whole_alignments_dup,
            "R1_chimeric_aligned_mates_dup" : self.r1_chimeric_alignments_dup,
            "R1_total_aligned_mates_dup" : self.r1_chimeric_alignments_dup +  self.r1_whole_alignments_dup,

            "R1_total_alignments_dedup" : self.r1_total_alignments_dedup,
            "R1_whole_aligned_mates_dedup" : self.r1_whole_alignments_dedup,
            "R1_chimeric_aligned_mates_dedup" : self.r1_chimeric_alignments_dedup,
            "R1_total_aligned_mates_dedup" : self.r1_chimeric_alignments_dedup +  self.r1_whole_alignments_dedup,

            "R2_total_alignments_dup" : self.r2_total_alignments_dup,
            "R2_whole_aligned_mates_dup" : self.r2_whole_alignments_dup,
            "R2_chimeric_aligned_mates_dup" : self.r2_chimeric_alignments_dup,
            "R2_total_aligned_mates_dup" : self.r2_chimeric_alignments_dup +  self.r2_whole_alignments_dup,
            
            "R2_total_alignments_dedup" : self.r2_total_alignments_dedup,
            "R2_whole_aligned_mates_dedup" : self.r2_whole_alignments_dedup,
            "R2_chimeric_aligned_mates_dedup" : self.r2_chimeric_alignments_dedup,
            "R2_total_aligned_mates_dedup" : self.r2_chimeric_alignments_dedup +  self.r2_whole_alignments_dedup,

            "discarded_alignments_overlap" : self.illegal_overlap_count,
            "discarded_alignments_post_trim" : self.illegal_post_trim_count,
            "pairs_with_multiple_valid_alignments" : self.at_least_two_alignments

        }

        stats_dict.update(self.contact_stats)
        
        stats_df = pd.DataFrame.from_dict(stats_dict, orient="index").T

        stats_df.to_csv(self.stats_path, index=False, sep="\t")
        
    def __init__(self, 
                 bam, 
                 out_prefix, 
                 chrom_sizes,
                 min_mapq=30, 
                 max_molecule_size=750, 
                 max_inter_align_gap=20):

        self.min_mapq = min_mapq
        self.max_molecule_size = max_molecule_size
        self.max_inter_align_gap = max_inter_align_gap

        self.contacts = f'{out_prefix}_contacts.txt.gz'
        self.contacts_dedup = f'{out_prefix}_contacts_dedup.txt.gz'
        self.contacts_short = f'{out_prefix}_contacts_dedup.short.gz'

        self.chimeras = f'{out_prefix}_chimeras.txt.gz'
        
        self.stats_path = f"{out_prefix}_contacts_trim_stats.txt"
                
        self.bam = bam
        self.masked_bam = f'{out_prefix}_trimmed.bam'

        self.process_bam()

        self.contact_stats = dedup_contacts(self.contacts, self.chimeras, self.contacts_dedup, save_raw = False)

        make_short(self.contacts_dedup, self.contacts_short, chrom_sizes)
        

        self.generate_stats()
