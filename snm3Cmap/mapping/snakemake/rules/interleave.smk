


def get_fastq(wildcards):
    out = {"interleaved" : f"{wildcards.id}_interleaved.fastq.gz"}
    return out