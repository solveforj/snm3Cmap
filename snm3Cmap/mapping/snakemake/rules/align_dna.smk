
if trim_output == "interleaved":

    rule align:
        input:
            rules.reformat.output
        output:
            temp("{id}.bam")
        threads: 
            config["align"]["threads"]
        params:
            reference_path=config["general"]["reference_path"],
            extra=config["align"]["bwa_interleaved_params"]
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return f"{wildcards.id}.bam"

elif trim_output == "separate":

    rule align_r1:
        input:
            rules.reformat_r1.output
        output:
            temp("{id}_R1.bam")
        threads: 
            config["align"]["threads"]
        params:
            reference_path=config["general"]["reference_path"],
            extra=config["align"]["bwa_R1_params"]
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """
    
    rule align_r2:
        input:
            rules.reformat_r2.output
        output:
            temp("{id}_R2.bam")
        threads: 
            config[modes]["align"]["threads"]
        params:
            reference_path=config["general"]["reference_path"],
            extra=config["align"]["bwa_R2_params"]
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return expand("{id}_R{mate}.bam", id=wildcards.id, mate=[1, 2])
