
mask_protocol = config["read_analysis"]["mask_protocol"]

if mask_protocol == "default":

    include: "mask.smk"

allc_protocol = config["read_analysis"]["allc_protocol"]

if allc_protocol == "default":

    if mode == "bsdna":

        include: "allc.smk"

    else:

        raise Exception("Allc cannot be generated for non-bisulfite-converted DNA.")



