
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
            extra=config["align"]["joint_params"]
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} -p {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return f"{wildcards.id}.bam"

elif trim_output == "separate" and not joint_alignments:

    rule align_r1:
        input:
            rules.reformat_r1.output
        output:
            temp("{id}_R1.bam")
        threads: 
            config["align"]["threads"]
        params:
            reference_path=config["general"]["reference_path"],
            extra=config["align"]["separate_R1_params"]
        retries: 3
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
            extra=config["align"]["separate_R2_params"]
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return expand("{id}_R{mate}.bam", id=wildcards.id, mate=[1, 2])

elif trim_output == "separate" and joint_alignments:

    rule align:
        input:
            get_trimmed_r1_fastq, 
            get_trimmed_r2_fastq
        output:
            temp("{id}.bam")
        threads: 
            config["align"]["threads"]
        params:
            reference_path=config["general"]["reference_path"],
            extra=config["align"]["joint_params"]
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return f"{wildcards.id}.bam"
