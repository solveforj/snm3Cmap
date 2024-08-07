
import subprocess
import shlex

rule get_dna_reads:
    input:
        unpack(get_fastq)
    output:
        dna_r1=temp("{id}_DNA_R1.fastq.gz"),
        dna_r2=temp("{id}_DNA_R2.fastq.gz"),
        rna_r1=temp("{id}_RNA_R1.fastq.gz"),
        rna_r2=temp("{id}_RNA_R2.fastq.gz"),
        stats=temp("{id}_cutadapt_stats.txt")
    threads:
        2
    params:
        cutadapt=config["trim_methods"]["hires"]["cutadapt"],
        extra=config["trim_methods"]["hires"]["cutadapt_params"]
    shell:
        """
        {params.cutadapt} --report=minimal {params.extra} -j {threads} \
            --untrimmed-output {output.dna_r1} --untrimmed-paired-output {output.dna_r2} \
            -o {output.rna_r1} -p {output.rna_r2} \
            {input.r1} {input.r2} > {output.stats}
        """

rule trim_stats:
    input:
        stats=rules.get_dna_reads.output.stats,
        r1=rules.get_dna_reads.output.dna_r1,
        r2=rules.get_dna_reads.output.dna_r2
    output:
        stats = temp("{id}_trim_stats.txt")
    threads:
        1
    params:
        seqtk=config["trim_methods"]["hires"]["seqtk"]
    run:
        with open(output["stats"], "w") as f, open(input["stats"]) as i:
            input_data = i.readlines()[1].strip().split()
            in_reads = input_data[1]
            out_rna_reads = input_data[6]
            
            seqtk = params["seqtk"]

            r1 = input["r1"]
            command1 = f"{seqtk} size {r1}"
            p1 = subprocess.run(shlex.split(command1), capture_output=True, text=True)
            r1_reads = p1.stdout.strip().split()[0]

            f.write("\t".join(["pre_trim_reads", "post_trim_rna_pairs", "post_trim_dna_pairs"]) + "\n")
            f.write("\t".join([in_reads, out_rna_reads, r1_reads]) + "\n")
    
trim_output = "separate"

def get_trimmed_r1_fastq(wildcards):
    return f"{wildcards.id}_DNA_R1.fastq.gz"

def get_trimmed_r2_fastq(wildcards):
    return f"{wildcards.id}_DNA_R2.fastq.gz"