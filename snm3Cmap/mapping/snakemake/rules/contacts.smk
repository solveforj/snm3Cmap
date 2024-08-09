

if config["contacts"]["call"]["call_protocol"] != "default":
    raise Exception("Contact calling protocol must be specified for desired output.")


rule generate_contacts:
    input:
        get_merged_bam
    output:
        bam="{id}_trimmed.bam",
        contacts=(
            "{id}_contacts.pairs.gz"
            if last_contacts_step == "call"
            else temp("{id}_contacts.pairs.gz")
        ),
        artefacts=(
            "{id}_artefacts.pairs.gz"
            if last_contacts_step == "call"
            else temp("{id}_artefacts.pairs.gz")
        ),
        stats=temp("{id}_alignment_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        chrom_sizes=config["general"]["chrom_sizes"],
        restriction_sites=config["contacts"]["call"]["restriction_sites"],
        extra=config["contacts"]["call"]["call_params"],
        manual_mate_annotation=('--manual-mate-annotation ' 
                                if trim_output == "separate" and not joint_alignments 
                                else ''),
        read_type="wgs" if mode == "dna" else "bisulfite"
        
    threads:
        1
    shell:
        'snm3Cmap call-contacts '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '--chrom-sizes {params.chrom_sizes} '
        '--restriction-sites {params.restriction_sites} '
        '--read-type {params.read_type} '
        '{params.manual_mate_annotation} '
        '{params.extra} '

def get_pairs_data(wildcards):
   return {"contacts" : f"{wildcards.id}_contacts.pairs.gz",
           "artefacts" : f"{wildcards.id}_artefacts.pairs.gz"}
    
def get_pairs_stats(wildcards):
   return {"contacts_stats": [],
           "artefacts_stats": []}

def get_filterbycov_stats(wildcards):
    return {"filterbycov_stats": []}

def get_trimmed_bam(wildcards):
    return f"{wildcards.id}_trimmed.bam"

if config["contacts"]["dedup"]["dedup_protocol"] == "default":

    include: "dedup_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_contacts.dedup.pairs.gz",
               "artefacts" : f"{wildcards.id}_artefacts.dedup.pairs.gz"
              }
        
    def get_pairs_stats(wildcards):
       return {"contacts_stats": f"{wildcards.id}_contacts_dedup_stats.txt",
               "artefacts_stats": f"{wildcards.id}_artefacts_dedup_stats.txt",
              }

if config["contacts"]["lowcov"]["lowcov_protocol"] == "default":

    include: "lowcov_contacts.smk"
    
    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_contacts.dedup.lowcov.pairs.gz",
               "artefacts" : f"{wildcards.id}_artefacts.dedup.highcov.pairs.gz"
              }
        
    def get_filterbycov_stats(wildcards):
        return {"filterbycov_stats": f"{wildcards.id}_filterbycov_stats.txt"}

rule pairtools_stats:
    input:
        unpack(get_pairs_data),
        unpack(get_pairs_stats),
        unpack(get_filterbycov_stats)
    output:
        temp("{id}_pairtools_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}"
    threads:
        1
    run:
        cmd = ("snm3Cmap pairtools-stats "
               "--out-prefix {params.out_prefix} "
               "--contacts {input.contacts} "
               "--artefacts {input.artefacts} "
              )
            
        if len(input["contacts_stats"]) != 0:
            cmd += "--contacts-stats {input.contacts_stats} "
        if len(input["artefacts_stats"]) != 0:
            cmd += "--artefacts-stats {input.artefacts_stats} "
        if len(input["filterbycov_stats"]) != 0:
            cmd += "--filterbycov-stats {input.filterbycov_stats} "

        shell(cmd)
        
rule coord_sort_trimmed:
    input:
        rules.generate_contacts.output.bam
    output:
        temp("{id}_trimmed_sorted.bam")
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

def get_coordsorted_analysis_bam(wildcards):
    return f"{wildcards.id}_trimmed_sorted.bam"
