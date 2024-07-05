import os

configfile: "config.yaml"

R1 = config["R1"]
R2 = config["R2"]
out_dir = config["out_dir"]
plate = config["plate"]
barcodes = config["barcodes"]

fq_out_r1 = f"{plate}_" "{cell_id}/" f"{plate}" "_{cell_id}_indexed_R1.fastq.gz"
fq_out_r2 = f"{plate}_" "{cell_id}/" f"{plate}" "_{cell_id}_indexed_R2.fastq.gz"

cutadapt_out_r1 = f"{plate}_" "{name}/" f"{plate}" "_{name}_indexed_R1.fastq.gz"
cutadapt_out_r2 = f"{plate}_" "{name}/" f"{plate}" "_{name}_indexed_R2.fastq.gz"

cell_ids = []
with open(barcodes) as f:
    for line in f:
        line = line.strip()
        if line[0] == ">":
            cell_ids.append(line[1:])

rule all:
    input:
        expand(fq_out_r1, cell_id=cell_ids),
        expand(fq_out_r2, cell_id=cell_ids),
        f"{plate}_demultiplex_stats.txt"

rule concatenate_R1:
    output:
        temp(f"{plate}_R1.fastq.gz")
    threads:
        1
    shell:
        'cat {R1} > {output}'

rule concatenate_R2:
    output:
        temp(f"{plate}_R2.fastq.gz")
    threads:
        1
    shell:
        'cat {R2} > {output}'


rule demultiplex:
    input:
        R1=rules.concatenate_R1.output,
        R2=rules.concatenate_R2.output
    output:
        R1=expand(fq_out_r1, cell_id=cell_ids),
        R2=expand(fq_out_r2, cell_id=cell_ids),
        stats=f"{plate}_demultiplex_stats.txt",
        R1_untrimmed="unknown_barcode_R1.fastq.gz",
        R2_untrimmed="unknown_barcode_R2.fastq.gz"
    threads:
        2
    shell:
        """
        cutadapt -j {threads} -Z -e 1 --no-indels --action=none -g ^file:{barcodes} \
            -o {cutadapt_out_r2}  \
            -p {cutadapt_out_r1} \
            --untrimmed-output {output.R2_untrimmed} \
            --untrimmed-paired-output {output.R1_untrimmed} \
            {input.R2} {input.R1} > {output.stats}
        """