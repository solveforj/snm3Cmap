from Bio.bgzf import BgzfWriter
from .read_trimmer import gap_pair_to_restriction_site
import re

class Pair:
    def __init__(self, algn1, 
                 algn2, readID,
                 ct, overlap, rule,
                 cs_locs, pair_index):
        
        self.algn1 = algn1
        self.algn2 = algn2
        self.readID = readID

        #print(ct)
        if "artefact" in ct:
            self.pair_class = "artefact"
        elif "na" != ct:
            self.pair_class = "contact"
        else:
            self.pair_class = ct

        
        self._rule = rule
        self._overlap = overlap
        self._ct = ct
        self._cs_locs = cs_locs
        self._pair_index = pair_index
        self._handle = None

        self.passed_filters = True
        
        self.chrom1 = f"\t{algn1['chrom']}"
        self.chrom2 = f"\t{algn2['chrom']}"
        self.pos1 = f"\t{str(algn1['pos'])}"
        self.pos2 = f"\t{str(algn2['pos'])}"
        self.strand1 = f"\t{algn1['strand']}"
        self.strand2 = f"\t{algn2['strand']}"
        
        self.phase0 = ""
        self.phase1 = ""
        self.pair_type = ""
        self.rule = ""
        self.reads = ""
        self.contact_class = ""
        self.multimap_overlap = ""
        self.cut_site_locs = ""

    def is_all(self):
        if self._rule == "all":
            return True
        elif self._rule == "mask":
            return False

    def add_phase(self, phase0, phase1):
        self.phase0 = f"\t{phase0}"
        self.phase1 = f"\t{phase1}"

    def add_metadata(self):
        
        self.pair_type = f"\t{self.algn1['type'] + self.algn2['type']}"
        self.rule = f"\t{self._rule}"
        self.reads = f"\t{self._pair_index[1]}"
        self.contact_class = f"\t{self._ct}"
        self.multimap_overlap = f"\t{str(self._overlap)}"
        self.cut_site_locs = f"\t{','.join(self._cs_locs)}"
        
    def __str__(self):
        line = (
            f"{self.readID}{self.chrom1}{self.pos1}{self.chrom2}{self.pos2}"
            f"{self.strand1}{self.strand2}{self.phase0}{self.phase1}"
            f"{self.pair_type}{self.rule}{self.reads}{self.contact_class}"
            f"{self.multimap_overlap}{self.cut_site_locs}\n"
        )
        return line
        
    
class PairsGenerator:
    def write_header(self, handle):
        
        handle.write("## pairs format v1.0\n")
        
        for chrom in self.chrom_sizes:
            handle.write(f'#chromsize: {chrom} {str(self.chrom_sizes[chrom])}\n')
        base_columns = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2"

        if self.snps:
            base_columns += " phase0 phase1"
        if self.full_pairs:
            base_columns += " pair_type rule reads contact_class multimap_overlap cut_site_locs"

        base_columns += "\n"
        handle.write(base_columns)

    def classify_pair(self, algn1,algn2, pair_index,
                         R1_trimmer,R2_trimmer,rule):

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

        contact_reads = pair_index[1]

        if rule == "all":
            # Pairtools reports 5' fragment before 3' fragment
            if contact_reads == "R1": 
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
            elif contact_reads == "R2": 
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
            elif contact_reads == "R1&2" or contact_reads == "R1-2":
                bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], 
                                                                                 algn2["read"], 
                                                                                 self.restriction_sites, 
                                                                                 self.max_cut_site_whole_algn_dist)
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
            bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], 
                                                                             algn2["read"], 
                                                                             self.restriction_sites, 
                                                                             self.max_cut_site_whole_algn_dist)
            if bp_enzyme != "artefact":
                ct = "gap"
            else:
                ct = "artefact_gap"
    
        return ct, overlap, cs_locs

    def pair_is_intra_short(self, pair):

        algn1 = pair.algn1
        algn2 = pair.algn2
        
        if algn1["chrom"] == algn2["chrom"]:
            if algn1["pos"] < algn2["pos"]:
                min_algn = algn1
                max_algn = algn2
            else:
                min_algn = algn2
                max_algn = algn1
    
            if max_algn["strand"] == min_algn["strand"]:
                if (max_algn["pos"] - min_algn["pos"]) < self.min_same_strand_dist:
                    pair.passed_filters = False
            elif min_algn["strand"] == "+":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_inward_dist:
                    pair.passed_filters = False
            elif min_algn["strand"] == "-":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_outward_dist:
                    pair.passed_filters = False

    def chrom_regex_pair(self, pair):

        if not self.chrom_regex.match(pair.algn1["chrom"]):
            pair.passed_filters = False

        if not self.chrom_regex.match(pair.algn2["chrom"]):
            pair.passed_filters = False
            
    def blacklist_pair(self, pair):
        
        pair.passed_filters = True

    def phase_pair(self, pair):

        pair.add_phase(".", ".")

    def metadata_pair(self, pair):

        pair.add_metadata()

    def all_pair(self, pair):

        if pair.is_all():
            pair.passed_filters = False

    def build_artefact_funcs(self):
        funcs = []
        funcs.append(self.pair_is_intra_short)
        
        if self.blacklist:
            funcs.append(self.blacklist_pair)
        if self.chrom_regex:
            funcs.append(self.chrom_regex_pair)
        #if self.snps:
        #    funcs.append(self.phase_pair)
        if self.full_pairs:
            funcs.append(self.metadata_pair)
        if self.remove_all:
            funcs.append(self.all_pair)
        self.artefact_funcs = funcs

    def build_contact_funcs(self):
        funcs = []
        funcs.append(self.pair_is_intra_short)
        
        if self.blacklist:
            funcs.append(self.blacklist_pair)
        if self.chrom_regex:
            funcs.append(self.chrom_regex_pair)
        if self.snps:
            funcs.append(self.phase_pair)
        if self.full_pairs:
            funcs.append(self.metadata_pair)
        if self.remove_all:
            funcs.append(self.all_pair)
        self.contact_funcs = funcs
            

    def write_pairs(self, algn1, algn2, readID, pair_index, R1_trimmer, R2_trimmer, rule):

        ct, overlap, cs_locs = self.classify_pair(algn1, algn2, pair_index, R1_trimmer, R2_trimmer, rule)

        #print(ct, overlap, cs_locs, readID)
        pair = Pair(algn1, algn2, readID, ct, overlap, rule, cs_locs, pair_index)
        #print(pair.pair_class)
        if pair.pair_class == "artefact":
                
            for filter in self.artefact_funcs:
                filter(pair)

            if pair.passed_filters:
                self.artefacts.write(str(pair))

        elif pair.pair_class == "contact":

            for filter in self.contact_funcs:
                filter(pair)
                #print(filter)
                #print(pair.passed_filters)

            if pair.passed_filters:
                self.contacts.write(str(pair))
        
        
    def __init__(self, 
                 contacts_path, 
                 chrom_sizes,
                 restriction_sites,
                 artefacts_path, 
                 blacklist=None,
                 snps=None,
                 remove_all=True,
                 full_pairs=False,
                 chrom_regex=None,
                 min_inward_dist=1000,
                 min_outward_dist=1000,
                 min_same_strand_dist=0,
                 max_cut_site_whole_algn_dist = 500
                ):
        
        self.artefacts_path = artefacts_path
        self.contacts_path = contacts_path
        self.chrom_sizes = chrom_sizes
        self.restriction_sites = restriction_sites
        self.blacklist = blacklist
        self.snps = snps
        self.remove_all = remove_all
        self.full_pairs = full_pairs
        self.chrom_regex = re.compile(chrom_regex) if chrom_regex else None
        self.min_inward_dist = min_inward_dist
        self.min_outward_dist = min_outward_dist
        self.min_same_strand_dist = min_same_strand_dist
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist

        self.build_contact_funcs()
        self.build_artefact_funcs()
        
    def __enter__(self):
        self.artefacts = BgzfWriter(self.artefacts_path, 'wb')
        self.write_header(self.artefacts)
            
        self.contacts = BgzfWriter(self.contacts_path, 'wb')
        self.write_header(self.contacts)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.artefacts.close()
        self.contacts.close()
        

