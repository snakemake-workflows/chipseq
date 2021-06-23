from smart_open import open
import os

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with open(snakemake.input[0]) as fin:
    with open(snakemake.output[0], 'w') as fout:
        for line in fin:
            items = line.split("\t")
            items[0] = os.path.basename(items[0])
            line = "{path}\t{col}".format(path=items[0], col=items[1])
            fout.write(line)
