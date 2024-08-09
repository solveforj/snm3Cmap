
rule dedup_contacts:
    input:
        contacts = rules.generate_contacts.output.contacts
    output:
        contacts_dedup = (
            "{id}_contacts.dedup.pairs.gz"
            if last_contacts_step == "dedup"
            else temp("{id}_contacts.dedup.pairs.gz")
        ),
        stats = temp("{id}_contacts_dedup_stats.txt")
    params:
        chrom_sizes=config["general"]["chrom_sizes"],
        flip_extra=config["contacts"]["dedup"]["flip_params"],
        sort_extra=config["contacts"]["dedup"]["sort_params"],
        dedup_extra=config["contacts"]["dedup"]["dedup_params"]
    threads:
        1
    shell:
        """
        contact_count=`zcat {input.contacts} | awk '!/^#/{{count++}} END{{ print count+0 }}'`
        if [ $contact_count -ge 1 ]; then
            pairtools flip {params.flip_extra} --chroms-path {params.chrom_sizes} {input.contacts} | \
            pairtools sort {params.sort_extra} --nproc {threads} - | \
            pairtools dedup {params.dedup_extra} -p {threads} --output {output.contacts_dedup} --output-stats {output.stats} - 
        else
            pairtools sort {params.sort_extra} --nproc {threads} --output {output.contacts_dedup} {input.contacts}
            rm {input.contacts}
            touch {output.stats}
        fi
        """

rule dedup_artefacts:
    input:
        artefacts = rules.generate_contacts.output.artefacts
    output:
        artefacts_dedup = (
            "{id}_artefacts.dedup.pairs.gz"
            if last_contacts_step == "dedup"
            else temp("{id}_artefacts.dedup.pairs.gz")
        ),
        stats = temp("{id}_artefacts_dedup_stats.txt")
    params:
        chrom_sizes=config["general"]["chrom_sizes"],
        flip_extra=config["contacts"]["dedup"]["flip_params"],
        sort_extra=config["contacts"]["dedup"]["sort_params"],
        dedup_extra=config["contacts"]["dedup"]["dedup_params"]
    threads:
        1
    shell:
        """
        contact_count=`zcat {input.artefacts} | awk '!/^#/{{count++}} END{{ print count+0 }}'`
        if [ $contact_count -ge 1 ]; then
            pairtools flip {params.flip_extra} --chroms-path {params.chrom_sizes} {input.artefacts} | \
            pairtools sort {params.sort_extra} --nproc {threads} - | \
            pairtools dedup {params.dedup_extra} -p {threads} --output {output.artefacts_dedup} --output-stats {output.stats} - 
        else
            pairtools sort {params.sort_extra} --nproc {threads} --output {output.artefacts_dedup} {input.artefacts} 
            rm {input.artefacts}
            touch {output.stats}
        fi
        """

