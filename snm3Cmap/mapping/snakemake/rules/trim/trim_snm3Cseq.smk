
rule trim_both:
    input:
        unpack(get_fastq)
    output:
        r1=temp("{id}_R1_trimmed.tmp.fastq.gz"),
        r2=temp("{id}_R2_trimmed.tmp.fastq.gz"),
        stats=temp("{id}_both_trim_stats.txt")
    threads: 
        1
    params:
        extra=config["preprocess"]["pair_trim_params"]
    shell:
        """
        cutadapt --report=minimal {params.extra} \
           -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} > {output.stats}
        """

rule trim_r1:
    input:
        rules.trim_both.output.r1
    output:
        fastq=temp("{id}_R1_trimmed.fastq.gz"),
        stats=temp("{id}_R1_trim_stats.txt")
    threads: 
        1
    params:
        extra=config["preprocess"]["r1_trim_params"]
    shell:
        """
        cutadapt --report=minimal {params.extra} \
            -o {output.fastq} {input} > {output.stats}
        """

rule trim_r2:
    input:
        rules.trim_both.output.r2
    output:
        fastq=temp("{id}_R2_trimmed.fastq.gz"),
        stats=temp("{id}_R2_trim_stats.txt")
    threads: 
        1
    params:
        extra=config["preprocess"]["r2_trim_params"]
    shell:
        """
        cutadapt --report=minimal {params.extra} \
            -o {output.fastq} {input} > {output.stats}
        """

rule merge_trim_stats:
    input:
        trim_both = rules.trim_both.output.stats,
        trim_r1 = rules.trim_r1.output.stats,
        trim_r2 = rules.trim_r2.output.stats
    output:
        stats = temp("{id}_trim_stats.txt")
    threads:
        1
    run:
        with open(input["trim_both"]) as f:
            lines = f.readlines()
            tstats = lines[1].strip().split("\t")
            in_pairs = tstats[1]
        with open(input["trim_r1"]) as f:
            lines = f.readlines()
            tstats = lines[1].strip().split("\t")
            out_r1 = tstats[6]
        with open(input["trim_r2"]) as f:
            lines = f.readlines()
            tstats = lines[1].strip().split("\t")
            out_r2 = tstats[6]
        with open(output["stats"], "w") as f:
            f.write("\t".join(["demultiplexed_pairs", "trimmed_R1_mates", "trimmed_R2_mates"]) + "\n")
            f.write("\t".join([in_pairs, out_r1, out_r2]) + "\n")

trim_output = "separate"

def get_trimmed_r1_fastq(wildcards):
    return f"{wildcards.id}_R1_trimmed.fastq.gz"

def get_trimmed_r2_fastq(wildcards):
    return f"{wildcards.id}_R2_trimmed.fastq.gz"
    