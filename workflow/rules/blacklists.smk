import os
import yaml
from smart_open import open

###### download igenomes file and blacklist ########

def remove_header(igenomes_link, igenomes_path):
    with open(igenomes_link) as fin:
        with open(igenomes_path,'w') as fout:
            for line in fin:
                if not line.strip().startswith('*') and not line.strip().startswith('/*') and not line.strip().startswith('//'):
                    fout.write(line)


def parse_to_yaml(igenomes):
    params = {"=": ":", " { ": ": {", "\"": "\'", " ": "", "fasta:": "\'fasta\':", "bwa:": "\'bwa\':",
              "bowtie2:": "\'bowtie2\':", "star:": "\'star\':", "bismark:": "\'bismark\':", "gtf:": "\'gtf\':",
              "bed12:": "\'bed12\':", "readme:": "\'readme\':", "mito_name:": "\'mito_name\':",
              "macs_gsize:": "\'macs_gsize\':", "blacklist:": "\'blacklist\':", "\'\'": "\', \'",
              "params:": "\'params\':", "genomes:": "\'genomes\':", ":": " : ", "{": " { ", "}": " } ",
              "} \'": "}, \'"}
    for i in params:
        igenomes = igenomes.replace(i,params[i])
    return igenomes


def add_links(igenomes):
    return igenomes.replace(
        "$ { baseDir } /","https://raw.githubusercontent.com/nf-core/chipseq/1.2.2/"
    ).replace(
        "$ { params.igenomes_base } /","s3://ngi-igenomes/igenomes/"
    )


def parse_igenomes(igenomes_link, igenomes_path):
    remove_header(igenomes_link,igenomes_path)
    with open(igenomes_path) as f:
        igenomes = yaml.load(add_links(parse_to_yaml(yaml.load(f,Loader=yaml.FullLoader))),Loader=yaml.FullLoader)
    with open(igenomes_path,'w') as f:
        yaml.dump(igenomes,f)


def generate_blacklist(build, chromosome, igenomes_path):
    with open(igenomes_path) as f:
        igenomes = yaml.load(f,Loader=yaml.FullLoader)
        if "blacklist" in igenomes["params"]["genomes"][build]:
            blacklist_link = igenomes["params"]["genomes"][build]["blacklist"]
            blacklist_path = get_blacklist_path(build, chromosome, igenomes_path, igenomes)
            with open(blacklist_link) as fin:
                with open(blacklist_path,'w') as fout:
                    for line in fin:
                        if chromosome:
                            if line.startswith("{}\t".format(chromosome)) or line.startswith("chr{}\t".format(chromosome)):
                                fout.write(line)
                        else:
                            fout.write(line)


def get_igenomes(igenomes_path):
    if os.path.isfile(igenomes_path):
        with open(igenomes_path) as f:
            return yaml.load(f,Loader=yaml.FullLoader)
    return None


def get_blacklist_path(build, chromosome, igenomes_path, igenomes):
    if config["resources"]["ref"]["blacklist"]:
        return config["resources"]["ref"]["blacklist"]
    elif "blacklist" in igenomes["params"]["genomes"][build]:
        if chromosome:
            return "{bl_dir}/chr{chr}_{build}-blacklist.bed".format(bl_dir=os.path.dirname(igenomes_path), chr=chromosome, build=build)
        return "{bl_dir}/{build}-blacklist.bed".format(bl_dir=os.path.dirname(igenomes_path), build=build)
    return ""
