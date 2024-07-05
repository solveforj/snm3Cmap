### 1. FASTQ files 

snm3Cseq profiles nuclei that are organized in 384-well plates, with one nuclei per well. It is required that each input FASTQ file contains data from only 1 plate. This is achieved in the Luo lab via demultiplexing using PCR indices with Illumina bcl2fastq. All FASTQ files for a plate should be in the same directory; however, this directory can contain FASTQ files for multiple plates. An example of the contents of this directory:

```
├── fastq/dir/
    ├── plateA_S1_L001_R1_001.fastq.gz
    ├── plateA_S1_L001_R2_001.fastq.gz
    ├── plateA_S1_L002_R1_001.fastq.gz
    ├── plateA_S1_L002_R2_001.fastq.gz
    ├── plateB_S1_L001_R1_001.fastq.gz
    ├── plateB_S1_L001_R2_001.fastq.gz
    ├── plateB_S1_L002_R1_001.fastq.gz
    └── plateB_S1_L002_R2_001.fastq.gz
```

### 2. plate_info.txt

This is a space/tab separated file with two columns. The first column contains the name of a plate. The second column contains the path to the FASTQ files to this plate's FASTQ files. For example:

```
plateA    /fastq/dir/
plateB    /fastq/dir/
```

### 3. barcodes

A FASTA file containing inline barcodes for snm3Cseq, which are the first 8 bp of R1. These determine what well a read pair came from in its plate. The latest snm3Cseq barcodes are provided in [random_index_v2.multiplex.fa](random_index_v2.multiplex.fa).

### 4. reference genome

The reference genome (i.e. hg38 or mm39) should be placed into a specific directory. Then, it should be indexed using Biscuit using the command below:

```{bash}
biscuit index /path/to/ref.fa
```

### 5. restriction cut sites

The location of the cut sites for the restriction enzymes used in the reference genome. snm3Cseq uses NlaIII and MboI, which will each have a file. These files can be generated using [this script](https://github.com/aidenlab/juicer/blob/main/misc/generate_site_positions.py) from Juicer. Specifically:
1. Modify the generate_site_positions.py script to include your reference genome path in the "filenames" variable
2. Modify the generate_site_positions.py script to include NlaIII and MboI cut sites in the "patterns" variable. MboI is already included. NlaIII's cut site is CATG.
3. Run the following commands:

```{bash}
python generate_site_positions.py MboI /path/to/ref.fa
python generate_site_positions.py NlaIII /path/to/ref.fa
```

This will generate two files:
1. ref_MboI.txt
2. ref_NlaIII.txt

### 5. chromsome sizes

A file showing chromosome sizes for the reference genome. An example can be found [here](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes)

### 7. mapping_config.yml

This is a YAML file that contains relevant parameters for both demultiplexing and mapping snm3Cseq. An documented example is provided in `mapping_config.yml`. This file will note the locations of all of the above input files. 
