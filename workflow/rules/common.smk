from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype = str).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

report: "../report/workflow.rst"

##### wildcard constraints #####

wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"])

####### helpers ###########

def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    fq2_present = pd.isnull(units.loc[(sample, unit), "fq2"])
    if isinstance(fq2_present, pd.core.series.Series):
        # if this is the case, get_fastqs cannot work properly
        raise ValueError(
            f"Multiple fq2 entries found for sample-unit combination {sample}-{unit}.\n"
            "This is most likely due to a faulty units.tsv file, e.g. "
            "a unit name is used twice for the same sample.\n"
            "Try checking your units.tsv for duplicates."
        )
    return fq2_present

def get_individual_fastq(wildcards):
    """Get individual raw FASTQ files from unit sheet, based on a read (end) wildcard"""
    if ( wildcards.read == "0" or wildcards.read == "1" ):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    elif wildcards.read == "2":
        return units.loc[ (wildcards.sample, wildcards.unit), "fq2" ]

def get_multiqc_input(wildcards):
    multiqc_input = []
    for (sample, unit) in units.index:
        reads = [ "1", "2" ]
        if is_single_end(sample, unit):
            reads = [ "0" ]
            multiqc_input.extend(expand (["logs/cutadapt/{sample}-{unit}.se.log"],
            sample = sample, unit = unit))
        else:
            multiqc_input.extend(expand (["logs/cutadapt/{sample}-{unit}.pe.log"],
            sample = sample, unit = unit))

        multiqc_input.extend(
            expand (
                [
                    "results/qc/fastqc/{sample}.{unit}.{reads}_fastqc.zip",
                    "results/qc/fastqc/{sample}.{unit}.{reads}.html",
                    "results/mapped/{sample}-{unit}.mapped.flagstat",
                    "results/mapped/{sample}-{unit}.mapped.idxstats",
                    "results/mapped/{sample}-{unit}.mapped.stats.txt"
                ],
                sample = sample,
                unit = unit,
                reads = reads
            )
        )
    for sample in samples.index:
        multiqc_input.extend(
            expand (
                [
                    "results/picard_dedup/{sample}.metrics.txt",
                    "results/picard_dedup/{sample}.picard_dedup.flagstat",
                    "results/picard_dedup/{sample}.picard_dedup.idxstats",
                    "results/picard_dedup/{sample}.picard_dedup.stats.txt",
                    "results/preseq/{sample}.lc_extrap",
                    "results/filtered/{sample}.filtered.flagstat",
                    "results/filtered/{sample}.filtered.idxstats",
                    "results/filtered/{sample}.filtered.stats.txt",
                    "results/orphan_rm_sorted/{sample}.orphan_rm_sorted.idxstats",
                    "results/orphan_rm_sorted/{sample}.orphan_rm_sorted.flagstat",
                    "results/orphan_rm_sorted/{sample}.orphan_rm_sorted.stats.txt",
                    "results/qc/multiple_metrics/{sample}.alignment_summary_metrics",
                    "results/qc/multiple_metrics/{sample}.base_distribution_by_cycle_metrics",
                    "results/qc/multiple_metrics/{sample}.base_distribution_by_cycle.pdf",
                    "results/qc/multiple_metrics/{sample}.insert_size_metrics",
                    "results/qc/multiple_metrics/{sample}.insert_size_histogram.pdf",
                    "results/qc/multiple_metrics/{sample}.quality_by_cycle_metrics",
                    "results/qc/multiple_metrics/{sample}.quality_by_cycle.pdf",
                    "results/qc/multiple_metrics/{sample}.quality_distribution_metrics",
                    "results/qc/multiple_metrics/{sample}.quality_distribution.pdf"
                ],
                sample = sample
            )
        )
    return multiqc_input

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def get_map_reads_input(wildcards):
    if is_single_end(wildcards.sample, wildcards.unit):
        return "results/trimmed/{sample}-{unit}.fastq.gz"
    return ["results/trimmed/{sample}-{unit}.1.fastq.gz", "results/trimmed/{sample}-{unit}.2.fastq.gz"]

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}-{unit}\tSM:{sample}-{unit}\tPL:{platform}'".format(
        sample=wildcards.sample,
        unit=wildcards.unit,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])
