
contamination_protocol = config["contamination"]["contamination_protocol"]

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
            extra=config["align"]["biscuit_interleaved_params"]
        retries: 3
        shell:
            """
            biscuit align -@ {threads} {params.extra} -M {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """

    if contamination_protocol == "default":
    
        rule bsconv:
            input:
                rules.align.output
            output:
                temp("{id}_bsconv.bam")
            threads:
                1
            params:
                reference_path=config["general"]["reference_path"]
            shell:
                """
                biscuit bsconv {params.reference_path} {input} {output} 
                """

        def get_bams(wildcards):
            return f"{wildcards.id}_bsconv.bam"
    
    else:
    
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
            extra=config["align"]["biscuit_R1_params"]
        retries: 3
        shell:
            """
            biscuit align -@ {threads} {params.extra} -M {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """
    
    rule align_r2:
        input:
            rules.reformat_r2.output
        output:
            temp("{id}_R2.bam")
        threads: 
            config["align"]["threads"]
        params:
            reference_path=config["general"]["reference_path"],
            extra=config["align"]["biscuit_R2_params"]
        retries: 3
        shell:
            """
            biscuit align -@ {threads} {params.extra} -M {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """
    
    
    if contamination_protocol == "default":
    
        rule bsconv_r1:
            input:
                rules.align_r1.output
            output:
                temp("{id}_R1_bsconv.bam")
            threads:
                1
            params:
                reference_path=config["general"]["reference_path"]
            shell:
                """
                biscuit bsconv {params.reference_path} {input} {output} 
                """
    
        rule bsconv_r2:
            input:
                rules.align_r2.output
            output:
                temp("{id}_R2_bsconv.bam")
            threads:
                1
            params:
                reference_path=config["general"]["reference_path"]
            shell:
                """
                biscuit bsconv {params.reference_path} {input} {output} 
                """
        
        def get_bams(wildcards):
            return expand("{id}_R{mate}_bsconv.bam", id=wildcards.id, mate=[1, 2])
    
    else:
    
        def get_bams(wildcards):
    
            return expand("{id}_R{mate}.bam", id=wildcards.id, mate=[1, 2])
