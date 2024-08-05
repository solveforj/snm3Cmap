if config["read_duplicates"]["duplicate_protocol"] == "default":

    rule coord_sort_mkdup:
        input:
            "{id}_mkdup.bam"
        output:
            temp("{id}_mkdup_sorted.bam")
        params:
            out_prefix=lambda wildcards: f"{wildcards.id}"
        threads: 
            10
        shell:
            """
            samtools sort -@ {threads} -o {output} {input}
            """

def get_all_stats(wildcards):
    stats_files = []
    if config["preprocess"]["trim_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_trim_stats.txt")
    if config["contamination"]["contamination_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_contam_stats.txt")
    if config["read_duplicates"]["duplicate_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_dupsifter_stats.txt")
    if config["align"]["align_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_alignment_stats.txt")
    if config["contacts"]["call_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_pairtools_stats.txt")  
    if config["read_analysis"]["allc_protocol"] != "none":
        stats_files.append(f"{wildcards.id}.allc.tsv.gz_methylation_stats.txt")  
    if config["read_duplicates"]["duplicate_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_mkdup_sorted.bam")  
    if config["contacts"]["call_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_trimmed_sorted.bam")  
    if config["read_analysis"]["mask_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_masked_sorted.bam")  

    return stats_files


rule aggregate_stats: 
    input:
        get_all_stats
    output:
        stats="{id}_qc_stats.txt"
    params:
        out_prefix = lambda wildcards: f"{wildcards.id}",
        mode = mode,
        extra = config["stats"]["stats_params"]
    threads:
        1
    shell:
        'snm3Cmap aggregate-qc-stats '
        '--job {params.out_prefix} '
        '--out-prefix {params.out_prefix} '
        '--mode {params.mode} '
        '{params.extra} '