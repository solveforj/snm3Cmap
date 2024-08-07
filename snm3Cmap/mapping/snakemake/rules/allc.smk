
rule bam_to_allc:
    input:
        get_coordsorted_analysis_bam
    output:
        allc = "{id}.allc.tsv.gz",
        tbi = "{id}.allc.tsv.gz.tbi",
        stats = "{id}.allc.tsv.gz.count.csv",
        methylation_stats = temp("{id}.allc.tsv.gz_methylation_stats.txt")
    params:
        reference_fasta=config["general"]["reference_path"],
        extra=config["read_analysis"]["allc"]["allc_params"]
    threads:
        1
    shell:
        'snm3Cmap bam-to-allc '
        '--bam-path {input} '
        '--reference-fasta {params.reference_fasta} '
        '--output-path {output.allc} ' 
        '--save-count-df ' 
        '{params.extra} '