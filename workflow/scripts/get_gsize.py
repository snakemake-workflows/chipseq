import yaml

def get_gsize_from_igenomes(igenomes, build):
    if build:
        with open(igenomes) as f:
            igenomes = yaml.load(f, Loader=yaml.FullLoader)
        if igenomes:
            if igenomes["params"]["genomes"][build]:
                if "macs_gsize" in igenomes["params"]["genomes"][build]:
                    return igenomes["params"]["genomes"][build]["macs_gsize"]
    return ""


igenomes = snakemake.input[0]
gsize_out = snakemake.output[0]

config_gsize = snakemake.params.get("extra", "")
build = snakemake.params.get("build", "")

if config_gsize:
    with open(gsize_out, 'w') as f:
        f.write("-g {}".format(config_gsize))
else:
    with open(gsize_out, 'w') as f:
        macs_gsize = get_gsize_from_igenomes(igenomes, build)
        if macs_gsize:
            f.write("-g {}".format(macs_gsize))
        else:
            f.write("")
