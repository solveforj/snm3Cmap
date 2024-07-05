# snm3C-seq

snm3C-seq jointly profiles methylation and chromatin conformation in single nuclei. The original method was described here:

Lee, Dong-Sung et al. “Simultaneous profiling of 3D genome structure and DNA methylation in single human cells.” Nature methods vol. 16,10 (2019): 999-1006. doi:10.1038/s41592-019-0547-z

## 1. Input

See [input_files.md](input_files.md) for a detailed description of the input files needed. 

## 2. Demultiplexing

Each well in 384-well plates contains a single nucleus. FASTQ files must be demultiplexed to produce FASTQ files that correspond to each nucleus.

### Run

First, run the following command

```{bash}
snm3Cmap prepare-demultiplex \
    --config /path/to/mapping_config.yml \
    --jobs 2 \
    --nolock \
    --rerun-incomplete
```

> For a more detailed description of the parameters of prepare-demultiplex, run
> ```{bash}
> snm3Cmap prepare-demultiplex --help
> ```

In the output directory specified in `mapping_config.yml` (`output/dir`), the following example structure will result:

```
├── output/dir/
    ├── plateA
        ├── plateA_A01/
        ├── plateA_A02/
        ├── ...
        ├── demultiplex_cmd.txt
        ├── config.yml
        └── demultiplex.smk
```
In the example above, the `/output/dir/plateA/plateA_A01/` directory will ultimately contain FASTQ files for nucleus in the A01 well of plateA. To perform demultiplexing of plateA, run the following command:

```{bash}
bash /output/dir/plateA/demultiplex_cmd.txt
```

Plates can be demultiplexed in parallel on a cluster.

### Output

The final output directory structure from the example above is shown below: 

```
├── output/dir/
    ├── plateA/
        ├── plateA_A01/
        │   ├── plateA_A01_indexed_R1.fastq.gz
        │   └── plateA_A01_indexed_R2.fastq.gz
        ├── plateA_A02/
        │   ├── plateA_A02_indexed_R1.fastq.gz
        │   └── plateA_A02_indexed_R2.fastq.gz
        ├── ...
        ├── plateA_demultiplex_stats.txt
        ├── unknown_barcode_R1.fastq.gz
        ├── unknown_barcode_R2.fastq.gz
        ├── demultiplex_cmd.txt
        ├── config.yml
        └── demultiplex.smk
```
In the example above, `plateA_A01_indexed_R1.fastq.gz` and `plateA_A01_indexed_R2.fastq.gz` represent R1/R2 FASTQ files for the nucleus in the A01 well of plateA.

Detailed documentation of how all output files can be interpreted can be found in [output_files.md](output_files.md).

## 2. Mapping

Requires demultiplexing to have been run, creating the output directory structure seen above.

### Run

First, run the following command:

```{bash}
snm3Cmap prepare-mapping \
    --config /path/to/mapping_config.yml \
    --jobs 30 \
    --nolock \
    --rerun-incomplete 
```

> For a more detailed description of the parameters of prepare-mapping, run
> ```{bash}
> snm3Cmap prepare-demultiplex --help
> ```
> 

In the output directory specified in `mapping_config.yml` (`output/dir`), the following example structure will result:

```
├── output/dir/
    ├── plateA/
        ├── plateA_A01/
        │   ├── plateA_A01_indexed_R1.fastq.gz
        │   ├── plateA_A01_indexed_R2.fastq.gz
        │   └── mapping_cmd.txt
        ├── plateA_A02/
        │   ├── plateA_A02_indexed_R1.fastq.gz
        │   ├── plateA_A02_indexed_R2.fastq.gz
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

To generate chromatin contacts and methylation data for wells A01 and A02 of plateA, run the following commands:

```{bash}
bash /output/dir/plateA/plateA_A01/mapping_cmd.txt
bash /output/dir/plateA/plateA_A02/mapping_cmd.txt
```

Well-level data can be mapped/processed in parallel on a cluster.

### Output

The final output directory structure from the example above is shown below: 

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
        ├── plateA_A02/
        │   ├── plateA_A02_indexed_R1.fastq.gz
        │   ├── plateA_A02_indexed_R2.fastq.gz
        │   ├── plateA_A02_artefacts.dedup.pairs.gz
        │   ├── plateA_A02_contacts.dedup.lowcov.pairs.gz
        │   ├── plateA_A02_trimmed.bam
        │   ├── plateA_A02.allc.tsv.gz
        │   ├── plateA_A02.allc.tsv.gz.count.csv
        │   ├── plateA_A02.allc.tsv.gz.tbi
        │   ├── plateA_A02_qc_stats.txt
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

In the example above, the contacts and methylation data for well A01 of plateA are `plateA_A01_contacts.dedup.lowcov.pairs.gz` and `plateA_A01.allc.tsv.gz`, respectively. Statistics that are relevant for quality control are found in `plateA_A01_qc_stats.txt`.

Detailed documentation of how all output files can be interpreted can be found in [output_files.md](output_files.md).

