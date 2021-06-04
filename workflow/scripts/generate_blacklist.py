import os
import yaml
from smart_open import open

# download blacklist and trim it for a specific chromosome


def copy_blacklist(igenomes, blacklist_path):
    with open(igenomes) as fin:
        with open(blacklist_path, 'w') as fout:
            for line in fin:
                fout.write(line)


def get_blacklist_from_igenomes(igenomes, blacklist_path):
    with open(igenomes) as f:
        igenomes = yaml.load(f, Loader=yaml.FullLoader)
        if "blacklist" in igenomes["params"]["genomes"][build]:
            blacklist_link = igenomes["params"]["genomes"][build]["blacklist"]
            with open(blacklist_link) as fin:
                with open(blacklist_path, 'w') as fout:
                    for line in fin:
                        if line.startswith("chr"):
                            line = line.replace("chr", "", 1)
                        if chromosome:
                            if line.startswith("{}\t".format(chromosome)):
                                fout.write(line)
                        else:
                            fout.write(line)
        else:
            open(blacklist_path, 'a').close()


igenomes = snakemake.input[0]
blacklist_path = snakemake.output.get("blacklist_path", "")

build = snakemake.params.get("build", "")
chromosome = snakemake.params.get("chromosome", "")
blacklist = snakemake.params.get("blacklist", "")

if blacklist:
    copy_blacklist(igenomes, blacklist_path)
else:
    get_blacklist_from_igenomes(igenomes, blacklist_path)
