
rule pairtools_filterbycov:
    input:
        unpack(get_pairs_data)
    output:
        lowcov = (
            "{id}_contacts.dedup.lowcov.pairs.gz"
            if last_contacts_step == "lowcov"
            else temp("{id}_contacts.dedup.lowcov.pairs.gz")
        ),
        highcov = temp("{id}_contacts.dedup.highcov.pairs.gz"),
        stats = temp("{id}_filterbycov_stats.txt"),
        merged = (
            "{id}_artefacts.dedup.highcov.pairs.gz"
            if last_contacts_step == "lowcov"
            else temp("{id}_artefacts.dedup.highcov.pairs.gz")
        )
    params:
        extra=config["contacts"]["lowcov"]["filterbycov_params"]
    threads:
        1
    shell:
        """
        contact_count=`zcat {input.contacts} | awk '!/^#/{{count++}} END{{ print count+0 }}' -`
        if [ $contact_count -ge 1 ]; then
            pairtools filterbycov \
                --output {output.lowcov} \
                --output-highcov {output.highcov} \
                --output-stats {output.stats} \
                {input.contacts}
            pairtools merge \
                --output {output.merged} \
                --nproc {threads} \
                {input.artefacts} {output.highcov}
        else
            cp {input.artefacts} {output.merged}
            touch {output.highcov}
            touch {output.stats}
            cp {input.contacts} {output.lowcov} 
        fi
        """

