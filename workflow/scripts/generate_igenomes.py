import os
import yaml
from smart_open import open

# download and parse igenomes file


def remove_header():
    with open(igenomes_link) as fin:
        with open(igenomes_path, 'w') as fout:
            for line in fin:
                if not line.strip().startswith('*'):
                    if not line.strip().startswith('/*'):
                        if not line.strip().startswith('//'):
                            fout.write(line)


def parse_to_yaml(igenomes):
    params = {"=": ":", " { ": ": {", "\"": "\'", " ": "", "fasta:": "\'fasta\':", "bwa:": "\'bwa\':",
              "bowtie2:": "\'bowtie2\':", "star:": "\'star\':", "bismark:": "\'bismark\':", "gtf:": "\'gtf\':",
              "bed12:": "\'bed12\':", "readme:": "\'readme\':", "mito_name:": "\'mito_name\':",
              "macs_gsize:": "\'macs_gsize\':", "blacklist:": "\'blacklist\':", "\'\'": "\', \'",
              "params:": "\'params\':", "genomes:": "\'genomes\':", ":": " : ", "{": " { ", "}": " } ",
              "} \'": "}, \'"}
    for i in params:
        igenomes = igenomes.replace(i, params[i])
    return igenomes


def add_links(igenomes):
    return igenomes.replace(
        "$ { baseDir } /", "https://raw.githubusercontent.com/nf-core/chipseq/1.2.2/"
    ).replace(
        "$ { params.igenomes_base } /", "s3://ngi-igenomes/igenomes/"
    )


def parse_igenomes():
    remove_header()
    with open(igenomes_path) as f:
        igenomes = yaml.load(add_links(parse_to_yaml(yaml.load(f, Loader=yaml.FullLoader))), Loader=yaml.FullLoader)
    with open(igenomes_path, 'w') as f:
        yaml.dump(igenomes, f)


igenomes_path = snakemake.output[0]
igenomes_release = snakemake.params.get("igenomes_release", "")
blacklist = snakemake.params.get("blacklist", "")

if igenomes_release:
    igenomes_link = "https://raw.githubusercontent.com/nf-core/chipseq/{version}/conf/igenomes.config".format(
        version=igenomes_release
    )
else:
    igenomes_link = "https://raw.githubusercontent.com/nf-core/chipseq/1.2.2/conf/igenomes.config"

parse_igenomes()
