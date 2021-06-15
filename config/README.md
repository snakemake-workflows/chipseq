
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the columns `sample`, `group`, `control`, and `antibody` have to be defined. 
* Samples / IP (immunoprecipitations) within the same `group` represents replicates and must have the same antibody and the same control.
* Controls / Input are listed like samples, but they do not have entries in the columns for `control` and `antibody`.
* The identifiers of each control has to be noted in the column `sample`.
* For all samples, the identifiers of the corresponding controls have to be given in the `control` column (see example below).

**Sample sheet example**:
* Samples / IP: A, B and C
* Controls / Input: D and E

| sample | group  | batch_effect | control | antibody |
|--------|--------|--------------|---------|----------|
| A      | TNFa   | batch1       | D       | p65      |
| B      | TNFa   | batch2       | D       | p65      |
| C      | E2TNFa | batch1       | E       | p65      |
| D      | TNFa   | batch1       |         |          |
| E      | E2TNFa | batch1       |         |          |

# Unit sheet

For each sample, add one or more sequencing units (runs or lanes) to the unit sheet `config/units.tsv`. For each unit, the columns `sample`, `unit`, `platform` and either `fq1` (single-end reads) or `fq1` and `fq2` (paired-end reads) or `sra_accession` have to be defined. 
* Each unit has a name specified in column `unit`, which can be e.g. a running number, or an actual run, lane or replicate id.
* Each unit has a `sample` name, which associates it with the biological sample it comes from.
* For single-end reads define for each unit either a path to FASTQ file (column `fq1`) or define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using the column `sra_accession`. 
* For paired-end reads define for each unit either two paths to FASTQ files (columns `fq1`, `fq2`) or define an SRA accession (column `sra_accession`).
* In case SRA accession numbers are used, the pipeline will automatically download the corresponding reads from SRA. If both local files and SRA accession are available, the local files will be preferred.
* The platform column needs to contain the used sequencing platform (one of 'CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'ONT', 'PACBIO’).

**Unit sheet example for single-end reads:**

| sample | unit | fq1                  | fq2 | sra_accession | platform |
|--------|------|----------------------|-----|---------------|----------|
| A      | 1    | data/A-run1.fastq.gz |     |               | ILLUMINA |
| B      | 1    | data/B-run1.fastq.gz |     |               | ILLUMINA |
| B      | 2    | data/B-run2.fastq.gz |     |               | ILLUMINA |
| C      | 1    | data/C-run1.fastq.gz |     |               | ILLUMINA |

**Unit sheet example for paired-end reads:**

| sample | unit | fq1                    | fq2                    | sra_accession | platform |
|--------|------|------------------------|------------------------|---------------|----------|
| A      | 1    | data/A-run1_1.fastq.gz | data/A-run1_2.fastq.gz |               | ILLUMINA |
| B      | 1    | data/B-run1_1.fastq.gz | data/B-run1_2.fastq.gz |               | ILLUMINA |
| B      | 2    | data/B-run2_1.fastq.gz | data/B-run1_2.fastq.gz |               | ILLUMINA |
| C      | 1    | data/C-run1_1.fastq.gz | data/C-run1_2.fastq.gz |               | ILLUMINA |

**Unit sheet example with SRA download (single-end or paired-end reads):**

| sample | unit | fq1 | fq2 | sra_accession | platform |
|--------|------|-----|-----|---------------|----------|
| A      | 1    |     |     | SRR1635456    | ILLUMINA |
| B      | 1    |     |     | SRR1635457    | ILLUMINA |
| B      | 2    |     |     | SRR1635458    | ILLUMINA |
| C      | 1    |     |     | SRR1635439    | ILLUMINA |
