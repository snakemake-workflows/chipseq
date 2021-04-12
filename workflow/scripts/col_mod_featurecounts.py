import sys
import pandas as pd
import os.path

sys.stderr = open(snakemake.log[0], "w")

fcounts = snakemake.input.get("featurecounts")
samples = pd.read_csv(snakemake.input.get("samples_file"), sep="\t")
bam = snakemake.input.get("bam", "")


def get_group_control_combination(bam_path):
    sample = os.path.basename(bam_path.split(".sorted.bam")[0])
    group = "".join(samples[samples["sample"] == sample]["group"])
    control = "".join(samples[samples["sample"] == sample]["control"])
    return "{}_{}.{}".format(group, control, sample)


def modify_header(old_header):
    return [get_group_control_combination(i) if i in bam else i for i in old_header]


f_counts_tab = pd.read_csv(fcounts, sep="\t", skiprows=1)
header = list(f_counts_tab.columns.values)
header_mod = modify_header(header)
f_counts_frame = pd.DataFrame(f_counts_tab.values)
f_counts_frame.columns = header_mod
f_counts_frame.to_csv(snakemake.output[0], sep='\t', index=False)

