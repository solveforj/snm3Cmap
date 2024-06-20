
import yaml

with open(parameters_path) as f:
    parameters = yaml.safe_load(f)

r1_fastq = config["r1"]
r2_fastq = config["r2"]
cell = config["cell"]

rule all:
    input:
        expand("{cell}_contacts.pairs", cell=cell),
        expand("{cell}_artefacts.pairs", cell=cell),
        expand("{cell}_trimmed.bam", cell=cell),
        expand("{cell}_alignment_stats.txt", cell=cell)
        
rule reformat:
    input:
        r1=r1_fastq,
        r2=r2_fastq
    output:
        r1=temp("{cell}_R1_trimmed_reformat.fastq.gz"),
        r2=temp("{cell}_R2_trimmed_reformat.fastq.gz")
    threads:
        1
    shell:
        """
        gzip -cd {input.r1} | awk -F' ' '{{l=l+1; if ((l-1)%4==0) {{$2=1;print $1"_"$2;}}else{{print $0}}}}' \
            | gzip -c > {output.r1} 
        
        gzip -cd {input.r2} | awk -F' ' '{{l=l+1; if ((l-1)%4==0) {{$2=2;print $1"_"$2;}}else{{print $0}}}}' \
            | gzip -c > {output.r2}
        """

rule align_r1:
    input:
        fastq=rules.reformat.output.r1
    output:
        bam_algn=temp("{cell}_R1.bam"),
    threads: 
        mapping_threads
    params:
        reference_path=parameters["align"]["reference_path"]
    shell:
        """
        bwa mem -5SPM -t {threads} {params.reference_path} {input.fastq} | samtools sort -o {output.bam_algn} -O BAM
        """

rule align_r2:
    input:
        fastq=rules.reformat.output.r2
    output:
        bam_algn=temp("{cell}_R2.bam"),
    threads: 
        mapping_threads
    params:
        reference_path=parameters["align"]["reference_path"]
    shell:
        """
        bwa mem -5SPM -t {threads} {params.reference_path} {input.fastq} | samtools sort -o {output.bam_algn} -O BAM
        """

rule merge_sort:
    input:
        r1=rules.align_r1.output.bam_algn,
        r2=rules.align_r2.output.bam_algn
    output:
        merged=temp("{cell}_merged.bam"),
        sorted=temp("{cell}_merged_sorted.bam")
    threads: 
        10
    shell:
        """
        samtools merge -@ {threads} -o {output.merged} {input.r1} {input.r2}
        samtools sort -@ {threads} -n -o {output.sorted} {output.merged}
        """

rule generate_contacts:
    input:
        bam = rules.merge_sort.output.sorted
    output:
        bam="{cell}_trimmed.bam",
        contacts="{cell}_contacts.pairs",
        artefacts="{cell}_artefacts.pairs",
        stats="{cell}_alignment_stats.txt"
    params:
        out_prefix=lambda wildcards: f"{wildcards.cell}",
        chrom_sizes=parameters["generate_contacts"]["chrom_sizes"],
        restriction_sites=parameters["generate_contacts"]["restriction_sites"],
        min_mapq=parameters["generate_contacts"]["min_mapq"],
        max_molecule_size=parameters["generate_contacts"]["max_molecule_size"],
        max_inter_align_gap=parameters["generate_contacts"]["max_inter_align_gap"],
        trim_reporting=parameters["generate_contacts"]["trim_reporting"],
        min_intra_distance=parameters["generate_contacts"]["min_intra_dist"],
        read_type=parameters["generate_contacts"]["read_type"],
        max_cut_site_split_algn_dist=parameters["generate_contacts"]["max_cut_site_split_algn_dist"],
        max_cut_site_whole_algn_dist=parameters["generate_contacts"]["max_cut_site_whole_algn_dist"]
    threads:
        1
    shell:
        'snm3Cmap call-contacts '
        '--bam {input.bam} '
        '--out-prefix {params.out_prefix} '
        '--chrom-sizes {params.chrom_sizes} '
        '--restriction-sites {params.restriction_sites} '
        '--min-mapq {params.min_mapq} '
        '--max-molecule-size {params.max_molecule_size} '
        '--max-inter-align-gap {params.max_inter_align_gap} '
        '--trim-reporting {params.trim_reporting} '
        '--min-intra-dist {params.min_intra_distance} '
        '--read-type {params.read_type} '
        '--max-cut-site-split-algn-dist {params.max_cut_site_split_algn_dist} '
        '--max-cut-site-whole-algn-dist {params.max_cut_site_whole_algn_dist} '
