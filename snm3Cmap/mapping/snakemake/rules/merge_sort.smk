
rule merge:
    input:
        get_bams
    output:
        temp("{id}_merged.bam")
    threads: 
        10
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} {input}
        """

rule sort:
    input:
        rules.merge.output
    output:
        temp("{id}_merged_sorted.bam")
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -n -o {output} {input}
        """

contamination_protocol = config["contamination"]["contamination_protocol"]

if contamination_protocol == "default":

    rule contam_filter:
        input:
            bam = rules.sort.output
        output:
            bam=temp("{id}_contam_filtered.bam"),
            stats=temp("{id}_contam_stats.txt")      
        params:
            out_prefix=lambda wildcards: f"{wildcards.id}",
            extra = config["contamination"]["params"]
        threads:
            1
        shell:
            'snm3Cmap contamination-filter --bam {input.bam} --out-prefix {params.out_prefix} {params.extra}'


    def get_merged_bam(wildcards):
        return f"{wildcards.id}_contam_filtered.bam"

else:

    def get_merged_bam(wildcards):
        return f"{wildcards.id}_merged_sorted.bam"