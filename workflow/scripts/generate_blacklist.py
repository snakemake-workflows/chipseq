import os
import yaml
from smart_open import open

# download blacklist and trim it for a specific chromosome

igenomes_path = snakemake.input[0]
blacklist_path = snakemake.output[0]

build = snakemake.params.get("build", "")
chromosome = snakemake.params.get("chromosome", "")

with open(igenomes_path) as f:
    igenomes = yaml.load(f, Loader=yaml.FullLoader)
    if "blacklist" in igenomes["params"]["genomes"][build]:
        blacklist_link = igenomes["params"]["genomes"][build]["blacklist"]
        blacklist_path = blacklist_path
        with open(blacklist_link) as fin:
            with open(blacklist_path, 'w') as fout:
                for line in fin:
                    if chromosome:
                        if line.startswith("{}\t".format(chromosome)):
                            fout.write(line)
                        elif line.startswith("chr{}\t".format(chromosome)):
                            fout.write(line)
                    else:
                        fout.write(line)
