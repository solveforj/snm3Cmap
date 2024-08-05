
rule interleave:
    input:
        unpack(get_fastq)
    output:
        temp("{id}_interleaved.fastq.gz")
    threads:
        1
    params:
        seqtk=config["preprocess"]["seqtk"]
    shell:
        """
        {params.seqtk} mergepe {input.r1} {input.r2} | gzip > {output}
        """

rule trim:
    input:
        rules.interleave.output
    output:
        temp("{id}_trimmed.fastq.gz")
    threads: 
        2
    params:
        pre_meta=config["preprocess"]["pre-meta"],
        extra=config["preprocess"]["pre-meta_params"]
    shell:
        """
        zcat {input} | {params.pre_meta} -t {threads} {params.extra} - | gzip > {output}
        """

rule trim_stats:
    input:
        rules.trim.output
    output:
        stats = temp("{id}_trim_stats.txt")
    params:
        seqtk=config["preprocess"]["seqtk"]
    threads:
        1
    run:
        with open(output["stats"], "w") as f:
            f.write("\t".join(["post_trim_reads", "post_trim_bp"]) + "\n")
        shell("{params.seqtk} size {input} >> {output.stats}")
            


trim_output = "interleaved"

def get_trimmed_interleaved_fastq(wildcards):
    return f"{wildcards.id}_trimmed.fastq.gz"