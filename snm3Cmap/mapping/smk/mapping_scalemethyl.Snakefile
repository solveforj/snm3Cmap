
import yaml

with open(parameters_path) as f:
    parameters = yaml.safe_load(f)

plate = config["plate"]
cell = config["cell"]

rule all:
    input:
        expand("{plate}_{cell}_contacts.dedup.lowcov.pairs.gz", plate=plate, cell=cell),
        expand("{plate}_{cell}_artefacts.dedup.pairs.gz", plate=plate, cell=cell),
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
        r2_temp2=temp("{plate}_{cell}_R2_trimmed2.tmp.fastq.gz"),
        stats=temp("{plate}_{cell}_trim_stats.txt")
    threads: 
        1
    params:
        r1_adapter=parameters["trim"]["r1_adapter"],
        r2_adapter=parameters["trim"]["r2_adapter"]
    shell:
        """
        cutadapt --report=minimal -u 10 -U 29 \
            -o {output.r1_temp} -p {output.r2_temp} \
            {input.r1} {input.r2} > {output.stats}

        cutadapt --report=minimal \
            -a {params.r1_adapter} -A {params.r2_adapter} \
            -o {output.r1} -p {output.r2_temp2} \
            {output.r1_temp} {output.r2_temp} >> {output.stats}

        cutadapt --report=minimal -u -10 \
            -o {output.r2} {output.r2_temp2} >> {output.stats}
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
        mapping_threads
    params:
        reference_path=parameters["align"]["reference_path"]
    shell:
        """
        biscuit align -@ {threads} -b 3 -M {params.reference_path} {input.fastq} \
            | samtools sort -o {output.bam_algn} -O BAM 

        biscuit bsconv {params.reference_path} {output.bam_algn} {output.bam_bsconv} 
        """

rule align_r2:
    input:
        fastq=rules.reformat.output.r2
    output:
        bam_algn=temp("{plate}_{cell}_R2.bam"),
        bam_bsconv=temp("{plate}_{cell}_R2_bsconv.bam")
    threads: 
        mapping_threads
    params:
        reference_path=parameters["align"]["reference_path"]
    shell:
        """
        biscuit align -@ {threads} -b 1 -M {params.reference_path} {input.fastq} \
            | samtools sort -o {output.bam_algn} -O BAM 

        biscuit bsconv {params.reference_path} {output.bam_algn} {output.bam_bsconv} 
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
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}",
        min_mapq=parameters["contam_filter"]["min_mapq"]
    threads:
        1
    shell:
        'snm3Cmap contamination-filter '
        '--bam {input.bam} '
        '--min-mapq {params.min_mapq} '
        '--out-prefix {plate}_{cell}'

rule mark_duplicates:
    input:
        bam = rules.contam_filter.output.bam
    output:
        bam = temp("{plate}_{cell}_mkdup.bam"),
        stats = temp("{plate}_{cell}_dupsifter_stats.txt")
    params:
        reference_path=parameters["mark_duplicates"]["reference_path"]
    threads:
        2
    shell:
        """
        dupsifter -s -o {output.bam} -O {output.stats} {params.reference_path} {input.bam} 
        """

rule generate_contacts:
    input:
        bam = rules.mark_duplicates.output.bam
    output:
        bam="{plate}_{cell}_trimmed.bam",
        contacts=temp("{plate}_{cell}_contacts.pairs"),
        artefacts=temp("{plate}_{cell}_artefacts.pairs"),
        stats=temp("{plate}_{cell}_alignment_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}",
        chrom_sizes=parameters["generate_contacts"]["chrom_sizes"],
        restriction_sites=parameters["generate_contacts"]["restriction_sites"],
        min_mapq=parameters["generate_contacts"]["min_mapq"],
        max_molecule_size=parameters["generate_contacts"]["max_molecule_size"],
        max_inter_align_gap=parameters["generate_contacts"]["max_inter_align_gap"],
        trim_reporting=parameters["generate_contacts"]["trim_reporting"],
        min_inward_distance=parameters["generate_contacts"]["min_inward_distance"],
        min_outward_distance=parameters["generate_contacts"]["min_outward_distance"],
        min_same_strand_distance=parameters["generate_contacts"]["min_same_strand_distance"],
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
        '--min-inward-dist {params.min_inward_distance} '
        '--min-outward-dist {params.min_outward_distance} '
        '--min-same-strand-dist {params.min_same_strand_distance} '
        '--read-type {params.read_type} '
        '--max-cut-site-split-algn-dist {params.max_cut_site_split_algn_dist} '
        '--max-cut-site-whole-algn-dist {params.max_cut_site_whole_algn_dist} '

rule pairtools_contacts:
    input:
        contacts = rules.generate_contacts.output.contacts
    output:
        contacts_dedup = temp("{plate}_{cell}_contacts.dedup.pairs.gz"),
        stats = temp("{plate}_{cell}_contacts_dedup_stats.txt")
    params:
        chrom_sizes=parameters["pairtools"]["chrom_sizes"],
        max_mismatch=parameters["pairtools"]["max_mismatch"]
    threads:
        1
    shell:
        """
        contact_count=`awk '!/^#/{{count++}} END{{ print count+0 }}' {input.contacts}`
        if [ $contact_count -ge 1 ]; then
            pairtools flip --chroms-path {params.chrom_sizes} {input.contacts} | \
            pairtools sort --nproc {threads} - | \
            pairtools dedup -p {threads} --max-mismatch {params.max_mismatch} \
                --output {output.contacts_dedup} --output-stats {output.stats} - 
        else
            bgzip -kf {input.contacts}
            pairtools sort --nproc {threads} --output {output.contacts_dedup} {input.contacts}.gz 
            rm {input.contacts}.gz
            touch {output.stats}
        fi
        """

rule pairtools_artefacts:
    input:
        artefacts = rules.generate_contacts.output.artefacts
    output:
        artefacts_dedup = temp("{plate}_{cell}_artefacts_temp.dedup.pairs.gz"),
        stats = temp("{plate}_{cell}_artefacts_dedup_stats.txt")
    params:
        chrom_sizes=parameters["pairtools"]["chrom_sizes"],
        max_mismatch=parameters["pairtools"]["max_mismatch"]
    threads:
        1
    shell:
        """
        contact_count=`awk '!/^#/{{count++}} END{{ print count+0 }}' {input.artefacts}`
        if [ $contact_count -ge 1 ]; then
            pairtools flip --chroms-path {params.chrom_sizes} {input.artefacts} | \
            pairtools sort --nproc {threads} - | \
            pairtools dedup -p {threads} --max-mismatch {params.max_mismatch} \
                --output {output.artefacts_dedup} --output-stats {output.stats} - 
        else
            bgzip -kf {input.artefacts}
            pairtools sort --nproc {threads} --output {output.artefacts_dedup} {input.artefacts}.gz 
            rm {input.artefacts}.gz
            touch {output.stats}
        fi
        """

rule pairtools_filterbycov:
    input:
        contacts_dedup = rules.pairtools_contacts.output.contacts_dedup,
        artefacts_dedup = rules.pairtools_artefacts.output.artefacts_dedup
    output:
        lowcov = "{plate}_{cell}_contacts.dedup.lowcov.pairs.gz",
        highcov = temp("{plate}_{cell}_contacts.dedup.highcov.pairs.gz"),
        stats = temp("{plate}_{cell}_filterbycov_stats.txt"),
        merged = "{plate}_{cell}_artefacts.dedup.pairs.gz"
    params:
        max_cov=parameters["pairtools_filterbycov"]["max_cov"],
        max_dist=parameters["pairtools_filterbycov"]["max_dist"],
        method=parameters["pairtools_filterbycov"]["method"]
    threads:
        1
    shell:
        """
        contact_count=`zcat {input.contacts_dedup} | awk '!/^#/{{count++}} END{{ print count+0 }}' -`
        if [ $contact_count -ge 1 ]; then
            pairtools filterbycov \
                --output {output.lowcov} \
                --output-highcov {output.highcov} \
                --output-stats {output.stats} \
                --max-cov {params.max_cov} \
                --max-dist {params.max_dist} \
                --method {params.method}  \
                {input.contacts_dedup}
            pairtools merge \
                --output {output.merged} \
                --nproc {threads} \
                {input.artefacts_dedup} {output.highcov}
        else
            cp {input.artefacts_dedup} {output.merged}
            touch {output.highcov}
            touch {output.stats}
            cp {input.contacts_dedup} {output.lowcov} 
        fi
        """

rule pairtools_stats:
    input:
        contacts = rules.pairtools_contacts.output.contacts_dedup,
        contacts_stats = rules.pairtools_contacts.output.stats,
        artefacts = rules.pairtools_artefacts.output.artefacts_dedup,
        artefacts_stats = rules.pairtools_artefacts.output.stats,
        filterbycov_stats = rules.pairtools_filterbycov.output.stats
    output:
        pairtools_stats = temp("{plate}_{cell}_pairtools_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}"
    threads:
        1
    shell:
        'snm3Cmap pairtools-stats '
        '--out-prefix {params.out_prefix} '
        '--contacts {input.contacts} '
        '--artefacts {input.artefacts} '
        '--contacts-stats {input.contacts_stats} '
        '--artefacts-stats {input.artefacts_stats} '
        '--filterbycov-stats {input.filterbycov_stats} '

rule mask:
    input:
        bam = rules.generate_contacts.output.bam
    output:
        bam = temp("{plate}_{cell}_masked.bam"),
    params:
        out_prefix=lambda wildcards: f"{wildcards.plate}_{wildcards.cell}",
        min_mapq=parameters["mask"]["min_mapq"]
    threads:
        1
    shell:
        'snm3Cmap mask-overlaps '
        '--bam {input.bam} '
        '--out-prefix {params.out_prefix} '
        '--min-mapq {params.min_mapq} '

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
    params:
        reference_fasta=parameters["bam_to_allc"]["reference_path"],
        num_upstr_bases=parameters["bam_to_allc"]["num_upstr_bases"],
        num_downstr_bases=parameters["bam_to_allc"]["num_downstr_bases"],
        min_mapq=parameters["bam_to_allc"]["min_mapq"],
        min_base_quality=parameters["bam_to_allc"]["min_base_quality"]
    threads:
        1
    shell:
        'snm3Cmap bam-to-allc '
        '--bam-path {input.bam} '
        '--reference-fasta {params.reference_fasta} '
        '--output-path {output.allc} ' 
        '--num-upstr-bases {params.num_upstr_bases} '
        '--num-downstr-bases {params.num_downstr_bases} '
        '--min-mapq {params.min_mapq} '
        '--min-base-quality {params.min_base_quality} '
        '--save-count-df ' 

rule aggregate_stats: 
    input:
        rules.trim.output.stats,
        rules.contam_filter.output.stats,
        rules.mark_duplicates.output.stats,
        rules.generate_contacts.output.stats,
        rules.pairtools_stats.output.pairtools_stats,
        rules.bam_to_allc.output.methylation_stats,
        rules.coord_sort_trimmed.output.bam,
        rules.coord_sort_mkdup.output.bam,
        rules.coord_sort_masked.output.bam
    output:
        stats="{plate}_{cell}_qc_stats.txt"
    params:
        out_prefix = lambda wildcards: f"{wildcards.plate}_{wildcards.cell}",
        min_mapq=parameters["aggregate_stats"]["min_mapq"],
        min_base_quality=parameters["aggregate_stats"]["min_base_quality"]
    threads:
        1
    shell:
        'snm3Cmap aggregate-qc-stats '
        '--cell {params.out_prefix} '
        '--out-prefix {params.out_prefix} '
        '--min-mapq {params.min_mapq} '
        '--min-base-quality {params.min_base_quality} '
