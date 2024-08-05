
rule mask:
    input:
        get_trimmed_bam
    output:
        temp("{id}_masked.bam"),
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        extra=config["read_analysis"]["mask_params"]
    threads:
        1
    shell:
        'snm3Cmap mask-overlaps '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '{params.extra} '

rule coord_sort_masked:
    input:
        rules.mask.output
    output:
        temp("{id}_masked_sorted.bam")
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

def get_coordsorted_analysis_bam(wildcards):
    return f"{wildcards.id}_masked_sorted.bam"