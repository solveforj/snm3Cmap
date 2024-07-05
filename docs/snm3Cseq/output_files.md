### Final output directory structure

After demultiplexing and mapping, the following directory structure should result for a single plate (plateA): 

```
├── output/dir/
    ├── plateA/
        ├── plateA_A01/
        │   ├── plateA_A01_indexed_R1.fastq.gz
        │   ├── plateA_A01_indexed_R2.fastq.gz
        │   ├── plateA_A01_artefacts.dedup.pairs.gz
        │   ├── plateA_A01_contacts.dedup.lowcov.pairs.gz
        │   ├── plateA_A01_trimmed.bam
        │   ├── plateA_A01.allc.tsv.gz
        │   ├── plateA_A01.allc.tsv.gz.count.csv
        │   ├── plateA_A01.allc.tsv.gz.tbi
        │   ├── plateA_A01_qc_stats.txt
        │   └── mapping_cmd.txt
        ├── ...
        ├── plateA_demultiplex_stats.txt
        ├── unknown_barcode_R1.fastq.gz
        ├── unknown_barcode_R2.fastq.gz
        ├── demultiplex_cmd.txt
        ├── config.yml
        ├── demultiplex.smk
        ├── mapping_scripts.txt
        └── mapping.smk 
```

### Plate-level output files

The following output files are generated for plateA:

* `plateA_demultiplex_stats.txt` - This contains demultiplexing stats from Cutadapt, whose [documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#cutadapt-s-output) is useful for interpreting it. Briefly, it reports how many reads were obtained for each barcode.
* `unknown_barcode_R1.fastq.gz` and `unknown_barcode_R2.fastq.gz` - These are read pairs that did not have a barcode detected in the first 8 bp of R1.
* `demultiplex_cmd.txt` - This is a script that can be run to demultiplex plateA.
* `config.yml` - Configuration file for Snakemake demultiplexing this plate.
* `demultiplex.smk` - Snakefile for demultiplexing this plate.
* `mapping_scripts.txt` - A txt file where each line is the path to a well's `mapping_cmd.txt`. This is useful for aggregating these paths for running mapping/processing in parallel.
* `mapping.smk` - A Snakefile that runs mapping/processing for any well in this plate on its demultiplexed FASTQ files.

### Well-level output files

* `plateA_A01_indexed_R1.fastq.gz` and `plateA_A01_indexed_R2.fastq.gz` - These are the well-level demultiplexed FASTQ files.
* `plateA_A01_artefacts.dedup.pairs.gz` - These are potential chromatin contacts that were ultimately classified as artefacts, either due to not being proximal to restriction sites or having [high contact coverage](https://pairtools.readthedocs.io/en/latest/cli_tools.html#pairtools-filterbycov). They have been deduplicated to remove PCR duplicates.
* `plateA_A01_contacts.dedup.lowcov.pairs.gz` - These are chromatin contacts that are proximal to restriction sites and do not fall into high contact coverage regions.
* `plateA_A01_trimmed.bam` - This is the BAM file containing alignments of the reads from Biscuit.
* `plateA_A01.allc.tsv.gz` - This is an [ALLC](https://lhqing.github.io/ALLCools/intro.html) file containing methylation data.
* `plateA_A01.allc.tsv.gz.count.csv` - This is a CSV that reports methylation counts for different methylation contexts.
* `plateA_A01.allc.tsv.gz.tbi` - This is a tabix index for the ALLC file.
* `plateA_A01_qc_stats.txt` - This is a TXT file that contains quality control statistics for the well. Full details of these statistics are described in the section below
* `mapping_cmd.txt` - THis is a script that can be run to map/process well A01 for plateA.

### QC statistics

There are many QC statistics, but the ones that are most relevant are defined below. Other QC stats are reported that are useful for understanding the performance of specific algorithms used in the pipeline. These shall be detailed, along with the algorithms, separately.

* cell - indicates the plate and well
* demultiplexed_pairs - Indicates the read pairs after demultiplexing
* trimmed_R1_mates - Indicates the number of R1 mates remaining after trimming
* trimmed_R2_mates - Indicates the number of R2 mates remaining after trimming
* R1_contam_pass - Indicates the number of R1 mates that passed the contamination check
* R1_contam_fail - Indicates the number of R1 mates that failed the contamination check
* R2_contam_pass - Indicates the number of R2 mates that passed the contamination check
* R2_contam_fail - Indicates the number of R2 mates that failed the contamination check
* pre_dedup_mates - Indicates the number of R1+R2 mates before deduplication with [dupsifter](https://github.com/huishenlab/dupsifter)
* duplicate_mates - Indicates the number of R1+R2 mates marked as duplicates by dupsifter
* alignment_dup_rate - Indicates the fraction of total R1+R2 mates that were marked as duplicates by dupsifter
* R1_total_alignments_dup - Indicates the number of R1 alignments (including duplicates)
* R1_whole_aligned_mates_dup - Indicates the number of R1 mates that were aligned without split alignments (including duplicates)
* R1_split_aligned_mates_dup - Indicates the number of R1 mates that were aligned with at least 1 split alignment (including duplicates)
* R1_total_aligned_mates_dup - Indicates the number of R1 mates that had at least 1 alignment, whether whole or split (including duplicates)
* R1_total_alignments_dedup - Indicates the number of R1 alignments (with duplicates)
* R1_whole_aligned_mates_dedup - Indicates the number of R1 mates that were aligned without split alignments (with duplicates)
* R1_split_aligned_mates_dedup - Indicates the number of R1 mates that were aligned with at least 1 split alignment (with duplicates)
* R1_total_aligned_mates_dedup - Indicates the number of R1 mates that had at least 1 alignment, whether whole or split (with duplicates)
* R2_total_alignments_dup - Indicates the number of R2 alignments (including duplicates)
* R2_whole_aligned_mates_dup - Indicates the number of R2 mates that were aligned without split alignments (including duplicates)
* R2_split_aligned_mates_dup - Indicates the number of R2 mates that were aligned with at least 1 split alignment (including duplicates)
* R2_total_aligned_mates_dup - Indicates the number of R2 mates that had at least 1 alignment, whether whole or split (including duplicates)
* R2_total_alignments_dedup - Indicates the number of R2 alignments (with duplicates)
* R2_whole_aligned_mates_dedup - Indicates the number of R2 mates that were aligned without split alignments (with duplicates)
* R2_split_aligned_mates_dedup - Indicates the number of R2 mates that were aligned with at least 1 split alignment (with duplicates)
* R2_total_aligned_mates_dedup - Indicates the number of R2 mates that had at least 1 alignment, whether whole or split (with duplicates)
* contacts_total - Total number of contacts
* contacts_dup_rate	- Duplicate rate of contacts from [pairtools dedup](https://pairtools.readthedocs.io/en/latest/cli_tools.html#pairtools-dedup)
* contacts_intra1kb	- Number of intrachromosomal contacts with a distance $ \geq $ 1kb
* contacts_intra10kb - Number of intrachromosomal contacts with a distance $ \geq $ 10kb
* contacts_intra20kb - Number of intrachromosomal contacts with a distance $ \geq $ 20kb
* contacts_inter - Number of interchromosomal contacts
* contacts_represented_read_pairs - Number of read pairs that had at least 1 contact called
* contacts_high_coverage_pairs - Number of contacts that were removed due to high coverage
* artefacts_total - Total number of artefacts
* artefacts_dup_rate - Duplicate rate of artefacts from pairtools dedup
* artefacts_intra1kb - Number of intrachromosomal artefacts with a distance $ \geq $ 1kb
* artefacts_intra10kb - Number of intrachromosomal artefacts with a distance $ \geq $ 10kb
* artefacts_intra20kb - Number of intrachromosomal artefacts with a distance $ \geq $ 20kb
* artefacts_inter - Number of interchromosomal artefacts with a distance $ \geq $ 20kb
* artefacts_represented_read_pairs - Number of read pairs that had at least 1 artefact called
* reference_coverage_dedup_mask - Number of bases in reference genome that were covered (excluding duplicates and overlapping portions of R1/R2)
* mapped_nucleotides_dedup_mask	- Number of mapped bases across reads (excluding duplicates and overlapping portions of R1/R2)
* mCG/CG_chrL - Percent of CG sites in chrL that are methylated
* mCH/CH_chrL - Percent of CH sites in chrL that are methylated
* mCG/CG_global - Percent of CG sites in non-chrL genome that are methylated
* mCH/CH_global - Percent of CH sites in non-chrL genome that are methylated
