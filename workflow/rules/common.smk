from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

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

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def is_control(sample):
    control = samples.loc[sample]["control"]
    return pd.isna(control) or pd.isnull(control)

def get_sample_control_combinations():
    sam_contr = []
    for sample in samples.index:
        if not is_control(sample):
            sam_contr.extend(expand(["{sample}-{control}"], sample = sample, control = samples.loc[sample]["control"]))
    return sam_contr

def get_sample_control_peak_combinations_list():
    sam_contr = []
    for sample in samples.index:
        if not is_control(sample):
            sam_contr.extend(expand(["{sample}-{control}.{peak}"], sample = sample, control = samples.loc[sample]["control"], peak = config["params"]["peak_analysis"]))
    return sam_contr

def get_peaks_count_plot_input():
    plot_in = []
    plot_in.extend(expand(["results/macs2_callpeak/peaks_count/{sam_contr_peak}.peaks_count.tsv"],
                          sam_contr_peak = get_sample_control_peak_combinations_list()))
    return plot_in

def get_frip_score_input():
    plot_in = []
    plot_in.extend(expand(["results/intersect/{sam_contr_peak}.peaks_frip.tsv"],
                          sam_contr_peak = get_sample_control_peak_combinations_list()))
    return plot_in

def get_plot_macs_qc_input():
    plot_in = []
    plot_in.extend(expand(["results/macs2_callpeak/{sam_contr_peak}_peaks.{peak}Peak"],
                          sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak_analysis"]))
    return plot_in

def get_plot_homer_annotatepeaks_input():
    plot_in = []
    plot_in.extend(expand(["results/homer/annotate_peaks/{sam_contr_peak}_peaks.annotatePeaks.txt"],
                          sam_contr_peak = get_sample_control_peak_combinations_list()))
    return plot_in

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
                    "results/qc/multiple_metrics/{sample}.quality_distribution.pdf",
                    "results/deeptools/plot_profile_data.tab",
                    "results/phantompeakqualtools/{sample}.phantompeak.spp.out",
                    "results/phantompeakqualtools/{sample}.spp_correlation_mqc.tsv",
                    "results/phantompeakqualtools/{sample}.spp_nsc_mqc.tsv",
                    "results/phantompeakqualtools/{sample}.spp_rsc_mqc.tsv"
                ],
                sample = sample
            )
        )
        if not is_control(sample):
            multiqc_input.extend(
                expand (
                    [
                        "results/deeptools/{sample}-{control}.fingerprint_qcmetrics.txt",
                        "results/deeptools/{sample}-{control}.fingerprint_counts.txt",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.xls",
                        "results/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.tsv",
                        "results/intersect/{sample}-{control}.{peak}.peaks_frip.tsv",
                        "results/homer/plots/mqc_{peak}_annotatepeaks.tsv"
                    ],
                sample = sample,
                control = samples.loc[sample]["control"],
                peak = config["params"]["peak_analysis"]
            )
        )
        if config["params"]["lc_extrap"]:
                multiqc_input.extend( expand(["results/preseq/{sample}.lc_extrap"], sample = sample))
    return multiqc_input

def all_input(wildcards):

    wanted_input = []

    # QC with fastQC and multiQC
    wanted_input.extend([
        "results/qc/multiqc/multiqc.html"
    ])

    # trimming reads
    for (sample, unit) in units.index:
        if is_single_end(sample, unit):
            wanted_input.extend(expand(
                    [
                        "results/trimmed/{sample}-{unit}.fastq.gz",
                        "results/trimmed/{sample}-{unit}.se.qc.txt"
                    ],
                    sample = sample,
                    unit = unit
                )
            )
        else:
            wanted_input.extend(
                expand (
                    [
                        "results/trimmed/{sample}-{unit}.1.fastq.gz",
                        "results/trimmed/{sample}-{unit}.2.fastq.gz",
                        "results/trimmed/{sample}-{unit}.pe.qc.txt"
                    ],
                    sample = sample,
                    unit = unit
            )
        )

    # mapping, merging and filtering bam-files
    for sample in samples.index:
        wanted_input.extend(
            expand (
                [
                    "results/IGV/big_wig/merged_library.bigWig.igv.txt",
                    "results/deeptools/plot_profile.pdf",
                    "results/deeptools/heatmap.pdf",
                    "results/deeptools/heatmap_matrix.tab",
                    "results/phantompeakqualtools/{sample}.phantompeak.pdf"
                ],
                sample = sample
            )
        )
        if not is_control(sample):
            wanted_input.extend(
                expand(
                    [
                        "results/deeptools/{sample}-{control}.plot_fingerprint.pdf",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_treat_pileup.bdg",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_control_lambda.bdg",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak",
                        "results/IGV/macs2_callpeak/{peak}/merged_library.{sample}-{control}.{peak}_peaks.igv.txt",
                        "results/homer/annotate_peaks/{sample}-{control}.{peak}_peaks.annotatePeaks.txt",
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_count.pdf",
                        # "results/intersect/plot_{peak}_peaks_frip_score.pdf",
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_macs2.pdf",
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_macs2_summary.txt",
                        "results/homer/plots/plot_{peak}_annotatepeaks.pdf"#,
                        # "results/homer/plots/plot_{peak}_annotatepeaks_summary.pdf"
                    ],
                sample = sample,
                control = samples.loc[sample]["control"],
                peak = config["params"]["peak_analysis"]
            )
        )
            if config["params"]["peak_analysis"] == "broad":
                wanted_input.extend(
                    expand(
                        [
                            "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.gappedPeak"
                        ],
                        sample = sample,
                        control = samples.loc[sample]["control"],
                        peak = config["params"]["peak_analysis"]
                    )
                )
            if config["params"]["peak_analysis"] == "narrow":
                wanted_input.extend(
                    expand(
                        [
                            "results/macs2_callpeak/{sample}-{control}.{peak}_summits.bed"

                        ],
                        sample = sample,
                        control = samples.loc[sample]["control"],
                        peak = config["params"]["peak_analysis"]
                    )
                )
    return wanted_input
