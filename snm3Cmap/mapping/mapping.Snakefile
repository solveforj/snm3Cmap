
plate = config["plate"]
cell = config["cell"]

r1_adapter = "AGATCGGAAGAGCACACGTCTGAAC"
r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGA"

rule all:
    input:
        expand("{plate}_{cell}_contacts.dedup.pairs.gz", plate=plate, cell=cell),
        expand("{plate}_{cell}_chimeras.dedup.pairs.gz", plate=plate, cell=cell),
        expand("{plate}_{cell}_contacts.short.gz", plate=plate, cell=cell),
        expand("{plate}_{cell}_trimmed.bam", plate=plate, cell=cell),
        expand("{plate}_{cell}.allc.tsv.gz", plate=plate, cell=cell),
        expand("{plate}_{cell}.allc.tsv.gz.tbi", plate=plate, cell=cell),
        expand("{plate}_{cell}.allc.tsv.gz.count.csv", plate=plate, cell=cell),
        expand("{plate}_{cell}_qc_stats.txt", plate=plate, cell=cell)

rule trim:
    input:
        r1="{plate}_{cell}_indexed_R1.fastq.gz",
        r2="{plate}_{cell}_indexed_R2.fastq.gz"
    output:
        r1=temp("{plate}_{cell}_R1_trimmed.fastq.gz"),
        r2=temp("{plate}_{cell}_R2_trimmed.fastq.gz"),
        r1_temp=temp("{plate}_{cell}_R1_trimmed.tmp.fastq.gz"),
        r2_temp=temp("{plate}_{cell}_R2_trimmed.tmp.fastq.gz"),
        stats=temp("{plate}_{cell}_trim_stats.txt")
    threads: 
        1
    shell:
        """
        cutadapt --report=minimal -q 20 -m 50 \
            -a {r1_adapter} -A {r2_adapter}  \
           -o {output.r1_temp} -p {output.r2_temp} \
            {input.r1} {input.r2} > {output.stats}
       
        cutadapt --report=minimal -u 18 -u -10 -m 30 \
            -o {output.r1} {output.r1_temp} >> {output.stats}
        
        cutadapt --report=minimal -u 10 -u -10 -m 30 \
            -o {output.r2} {output.r2_temp} >> {output.stats}
        """
        
rule reformat:
    input:
        r1=rules.trim.output.r1,
        r2=rules.trim.output.r2
    output:
        r1=temp("{plate}_{cell}_R1_trimmed_reformat.fastq.gz"),
        r2=temp("{plate}_{cell}_R2_trimmed_reformat.fastq.gz")
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
        bam_algn=temp("{plate}_{cell}_R1.bam"),
        bam_bsconv=temp("{plate}_{cell}_R1_bsconv.bam")
    threads: 
        biscuit_threads
    shell:
        """
        biscuit align -@ {threads} -b 3 -M {reference_path} {input.fastq} \
            | samtools sort -o {output.bam_algn} -O BAM 

        biscuit bsconv {reference_path} {output.bam_algn} {output.bam_bsconv} 
        """

rule align_r2:
    input:
        fastq=rules.reformat.output.r2
    output:
        bam_algn=temp("{plate}_{cell}_R2.bam"),
        bam_bsconv=temp("{plate}_{cell}_R2_bsconv.bam")
    threads: 
        biscuit_threads
    shell:
        """
        biscuit align -@ {threads} -b 1 -M {reference_path} {input.fastq} \
            | samtools sort -o {output.bam_algn} -O BAM 

        biscuit bsconv {reference_path} {output.bam_algn} {output.bam_bsconv} 
        """

rule merge_sort:
    input:
        r1=rules.align_r1.output.bam_bsconv,
        r2=rules.align_r2.output.bam_bsconv
    output:
        merged=temp("{plate}_{cell}_merged.bam"),
        sorted=temp("{plate}_{cell}_merged_sorted.bam")
    threads: 
        10
    shell:
        """
        samtools merge -@ {threads} -o {output.merged} {input.r1} {input.r2}
        samtools sort -@ {threads} -n -o {output.sorted} {output.merged}
        """

rule contam_filter:
    input:
        bam = rules.merge_sort.output.sorted
    output:
        bam=temp("{plate}_{cell}_contam_filtered.bam"),
        stats=temp("{plate}_{cell}_contam_stats.txt")      
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads:
        1
    shell:
        'snm3Cmap contamination-filter '
        '--bam {input.bam} '
        '--mapq-threshold 30 '
        '--out-prefix {plate}_{cell}'

rule mark_duplicates:
    input:
        bam = rules.contam_filter.output.bam
    output:
        bam = temp("{plate}_{cell}_mkdup.bam"),
        stats = temp("{plate}_{cell}_dupsifter_stats.txt")
    threads:
        2
    shell:
        """
        dupsifter -s -o {output.bam} -O {output.stats} {reference_path} {input.bam} 
        """

rule generate_contacts:
    input:
        bam = rules.mark_duplicates.output.bam
    output:
        bam="{plate}_{cell}_trimmed.bam",
        contacts="{plate}_{cell}_contacts.dedup.pairs.gz",
        chimeras="{plate}_{cell}_chimeras.dedup.pairs.gz",
        short="{plate}_{cell}_contacts.short.gz",
        stats=temp("{plate}_{cell}_alignment_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads:
        1
    shell:
        'snm3Cmap call-contacts '
        '--bam {input.bam} '
        '--out-prefix {params.out_prefix} '
        '--chrom-sizes {chrom_sizes} '
        '--min-mapq {min_mapq} '
        '--max-molecule-size {max_molecule_size} '
        '--max-inter-align-gap {max_inter_align_gap} '
        '{full_bam} '
        '--min-intra-dist {min_intra_dist} '
        '--threads {threads} '

rule mask:
    input:
        bam = rules.generate_contacts.output.bam
    output:
        bam = temp("{plate}_{cell}_masked.bam"),
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads:
        1
    shell:
        'snm3Cmap mask-overlaps '
        '--bam {input.bam} '
        '--out-prefix {params.out_prefix} '
        '--min-mapq {min_mapq} '

rule coord_sort_trimmed:
    input:
        bam = rules.generate_contacts.output.bam
    output:
        bam = temp("{plate}_{cell}_trimmed_sorted.bam")
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        """

rule coord_sort_mkdup:
    input:
        bam = rules.mark_duplicates.output.bam
    output:
        bam = temp("{plate}_{cell}_mkdup_sorted.bam")
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        """

rule coord_sort_masked:
    input:
        bam = rules.mask.output.bam
    output:
        bam = temp("{plate}_{cell}_masked_sorted.bam")
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        """

rule bam_to_allc:
    input:
        bam = rules.coord_sort_masked.output.bam
    output:
        allc = "{plate}_{cell}.allc.tsv.gz",
        tbi = "{plate}_{cell}.allc.tsv.gz.tbi",
        stats = "{plate}_{cell}.allc.tsv.gz.count.csv",
        methylation_stats = temp("{plate}_{cell}.allc.tsv.gz_methylation_stats.txt")
    threads:
        1
    shell:
        'snm3Cmap bam-to-allc '
        '--bam-path {input.bam} '
        '--reference-fasta {reference_path} '
        '--output-path {output.allc} ' 
        '--min-mapq {min_mapq} '
        '--min-base-quality {min_base_quality} '
        '--save-count-df ' 

rule aggregate_stats: 
    input:
        rules.trim.output.stats,
        rules.contam_filter.output.stats,
        rules.mark_duplicates.output.stats,
        rules.generate_contacts.output.stats,
        rules.bam_to_allc.output.methylation_stats,
        rules.coord_sort_trimmed.output.bam,
        rules.coord_sort_mkdup.output.bam,
        rules.coord_sort_masked.output.bam
    output:
        stats="{plate}_{cell}_qc_stats.txt"
    params:
        out_prefix = lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads:
        1
    shell:
        'snm3Cmap aggregate-qc-stats '
        '--cell {params.out_prefix} '
        '--out-prefix {params.out_prefix} '
        '--min-mapq {min_mapq} '
        '--min-base-quality {min_base_quality} '
