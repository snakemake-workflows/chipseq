import sys
import pandas as pd
import os.path

sys.stderr = open(snakemake.log[0], "w")

fcounts = snakemake.input.get("featurecounts")
bam = snakemake.input.get("bam", "")

# common_prefix = os.path.commonprefix(bam)
# common_suffix = os.path.commonprefix([i[::-1] for i in bam])[::-1]

def strip_path(path):
    return os.path.splitext(os.path.basename(path))[0]


def modify_header(header):
    return [strip_path(i) if i in bam else i for i in header]

f_counts_tab = pd.read_csv(fcounts, sep="\t", skiprows=1)
header = list(f_counts_tab.columns.values)
header_mod = modify_header(header)
f_counts_frame = pd.DataFrame(f_counts_tab.values)
f_counts_frame.columns = header_mod
f_counts_frame.to_csv(snakemake.output[0], sep='\t', index=False)
