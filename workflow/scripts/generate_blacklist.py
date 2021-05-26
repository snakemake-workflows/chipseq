import os
import re
import yaml
from smart_open import open

# download blacklist and trim it for a specific chromosome


def copy_blacklist(igenomes_or_blacklist, blacklist_path):
    with open(igenomes_or_blacklist) as fin:
        with open(blacklist_path, 'w') as fout:
            for line in fin:
                fout.write(line)


def get_blacklist_from_igenomes(igenomes_or_blacklist, blacklist_path):
    with open(igenomes_or_blacklist) as f:
        igenomes = yaml.load(f, Loader=yaml.FullLoader)
        if "blacklist" in igenomes["params"]["genomes"][build]:
            blacklist_link = igenomes["params"]["genomes"][build]["blacklist"]
            with open(blacklist_link) as fin:
                with open(blacklist_path, 'w') as fout:
                    for line in fin:
                        if chromosome:
                            if line.startswith("{}\t".format(chromosome)):
                                fout.write(line)
                            elif line.startswith("chr{}\t".format(chromosome)):
                                pattern = re.compile(r'^chr')
                                fout.write(re.sub(pattern, "", line))
                        else:
                            pattern = re.compile(r'^chr')
                            fout.write(re.sub(pattern, "", line))
        else:
            open(blacklist_path, 'a').close()


igenomes_or_blacklist = snakemake.input[0]
blacklist_path = snakemake.output.get("blacklist_path", "")

build = snakemake.params.get("build", "")
chromosome = snakemake.params.get("chromosome", "")
blacklist = snakemake.params.get("blacklist", "")

if blacklist:
    copy_blacklist(igenomes_or_blacklist, blacklist_path)
else:
    get_blacklist_from_igenomes(igenomes_or_blacklist, blacklist_path)




