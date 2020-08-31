rule preseq_lc_extrap:
    input:
        "results/sam-view/{sample}.bam"
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        "-v -seed 1"   #optional parameters
    log:
        "logs/preseq/{sample}.log"
    wrapper:
        "0.64.0/bio/preseq/lc_extrap"

rule collect_multiple_metrics:
    input:
         bam="results/orphan_rm_sorted/{sample}.bam",
         ref="resources/ref/genome.fasta"
    output:
        # Through the output file extensions the different tools for the metrics can be selected
        # so that it is not necessary to specify them under params with the "PROGRAM" option.
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html.
        multiext("{path}{sample}",
                 ".alignment_summary_metrics",
                 ".base_distribution_by_cycle_metrics",
                 ".base_distribution_by_cycle.pdf",
                 ".insert_size_metrics",
                 ".insert_size_histogram.pdf",
                 ".quality_by_cycle_metrics",
                 ".quality_by_cycle.pdf",
                 ".quality_distribution_metrics",
                 ".quality_distribution.pdf"
                 )
    resources:
        # This parameter (default 3 GB) can be used to limit the total resources a pipeline is allowed to use, see:
        #     https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources
        mem_gb=3
    log:
        "logs/picard/{path}{sample}.log"
    params:
        # optional parameters
        "VALIDATION_STRINGENCY=LENIENT "
    wrapper:
        "0.64.0/bio/picard/collectmultiplemetrics"

rule genomecov:
    input:
        "results/orphan_rm_sorted/{sample}.bam",
        "results/orphan_rm_sorted/{sample}.orphan_rm_sorted.flagstat"
    output:
        pipe("results/bed_graph/{sample}.bedgraph")
    log:
        "logs/bed_graph/{sample}.log"
    params: #-fs option was not used because there are no single end reads any more
        "-bg -pc -scale $(grep 'mapped (' results/orphan_rm_sorted/{sample}.orphan_rm_sorted.flagstat | awk '{print 1000000/$1}')"
    wrapper:
        "0.64.0/bio/bedtools/genomecov"

rule sort_genomecov:
    input:
        "results/bed_graph/{sample}.bedgraph"
    output:
        "results/bed_graph/{sample}.sorted.bedgraph"
    log:
        "logs/sort_genomecov/{sample}.log"
    shell:
        "sort -k1,1 -k2,2n {input} > {output} 2> {log}"

rule bedGraphToBigWig:
    input:
        bedGraph="results/bed_graph/{sample}.sorted.bedgraph",
        chromsizes="resources/ref/genome.chrom.sizes"
    output:
        "results/big_wig/{sample}.bigWig"
    log:
        "logs/big_wig/{sample}.log"
    params:
        ""
    wrapper:
        "0.64.0/bio/ucsc/bedGraphToBigWig"

rule create_igv:
    input:
        "resources/ref/genome.bed",
        expand("results/big_wig/{sample}.bigWig", sample=samples.index)
    output:
        "results/IGV/merged_library.bigWig.igv.txt"
    log:
        "logs/igv/create_igv.log"
    shell:
        "find {input} -type f -name '*.bigWig' -exec echo -e 'results/big_wig/\"{{}}\"\t0,0,178' \;  > {output} 2> {log}"

rule compute_matrix:
    input:
         bed="resources/ref/genome.bed",
         bigwig=expand("results/big_wig/{sample}.bigWig", sample=samples.index)
    output:
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz="results/deeptools/matrix_files/matrix.gz",
        matrix_tab="results/deeptools/matrix_files/matrix.tab"
    log:
        "logs/deeptools/compute_matrix.log"
    params:
        command="scale-regions",
        extra="--regionBodyLength 1000 "
              "--beforeRegionStartLength 3000 "
              "--afterRegionStartLength 3000 "
              "--skipZeros "
              "--smartLabels "
              "--numberOfProcessors 2 "
    wrapper:
        "0.64.0/bio/deeptools/computematrix"

rule plot_profile:
    input:
         "results/deeptools/matrix_files/matrix.gz"
    output:
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotprofile.html.
        plot_img="results/deeptools/plot_profile.pdf",
        data="results/deeptools/plot_profile_data.tab"
    log:
        "logs/deeptools/plot_profile.log"
    params:
        ""
    wrapper:
        "0.64.0/bio/deeptools/plotprofile"

rule plot_heatmap:
    input:
         "results/deeptools/matrix_files/matrix.gz"
    output:
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotheatmap.html.
        heatmap_img="results/deeptools/heatmap.pdf",
        heatmap_matrix="results/deeptools/heatmap_matrix.tab"
    log:
        "logs/deeptools/heatmap.log"
    params:
        ""
    wrapper:
        "0.64.0/bio/deeptools/plotheatmap"

rule phantompeakqualtools:
    input:
         "results/orphan_rm_sorted/{sample}.bam"
    output:
        res_phantom="results/phantompeakqualtools/{sample}.phantompeak.spp.out",
        r_data="results/phantompeakqualtools/{sample}.phantompeak.Rdata",
        plot="results/phantompeakqualtools/{sample}.phantompeak.pdf"
    threads:
        8
    log:
        "logs/phantompeakqualtools/{sample}.phantompeak.log"
    conda:
        "../envs/phantompeakqualtools.yaml"
    shell:
        "( Rscript -e \"library(caTools); source('../workflow/scripts/run_spp.R')\" "
        "  -c={input} -savp={output.plot} -savd={output.r_data} "
        "  -out={output.res_phantom} -p={threads} 2>&1 ) >{log}"

rule phantompeak_correlation:
    input:
        data="results/phantompeakqualtools/{sample}.phantompeak.Rdata",
        header="../workflow/header/spp_corr_header.txt"
    output:
        "results/phantompeakqualtools/{sample}.spp_correlation_mqc.tsv"
    log:
        "logs/phantompeakqualtools/correlation/{sample}.spp_corr.log"
    script:
        "../scripts/phantompeak_correlation.R"

rule phantompeak_multiqc:
    # NSC (Normalized strand cross-correlation) and RSC (relative strand cross-correlation) metrics use cross-correlation
    input:
        data="results/phantompeakqualtools/{sample}.phantompeak.spp.out",
        nsc_header="../workflow/header/nsc_header.txt",
        rsc_header="../workflow/header/rsc_header.txt"
    output:
        nsc="results/phantompeakqualtools/{sample}.spp_nsc_mqc.tsv",
        rsc="results/phantompeakqualtools/{sample}.spp_rsc_mqc.tsv"
    log:
        "logs/phantompeakqualtools/correlation/{sample}.nsc_rsc.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "( gawk -v OFS='\t' '{{print $1, $9}}' {input.data} | cat {input.nsc_header} - > {output.nsc} && "
        "  gawk -v OFS='\t' '{{print $1, $10}}' {input.data} | cat {input.rsc_header} - > {output.rsc} 2>&1 ) >{log}"
