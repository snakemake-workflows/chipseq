from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

out_counts = snakemake.output.get("counts")

optional_output = ""

if out_counts:
    optional_output += " --outRawCounts {out_counts} ".format(out_counts=out_counts)

shell(
    "(plotFingerprint "
    "-b {snakemake.input.bam_files} "
    "-o {snakemake.output.fingerprint} "
    "{optional_output} "
    "{snakemake.params}) {log}"
)
