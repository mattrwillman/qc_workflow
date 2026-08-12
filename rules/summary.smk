# Cross-sample summaries.

rule coverage_summary:
    input:
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.coverage_stats.tsv"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
    output:
        "results/summary/coverage_summary.tsv",
    shell:
        """
        mkdir -p $(dirname {output})
        (head -n1 {input[0]}; for f in {input}; do tail -n +2 "$f"; done) > {output}
        """

rule snp_density_heatmap:
    # All samples, one heatmap per (aligner, param_set): rows=sample,
    # columns=genomic window, color=SNP count.
    input:
        snp_density_tsvs_for
    output:
        "results/mapping/{aligner}/{param_set}/variants/snp_density_heatmap.png",
    log:
        "logs/summary/{aligner}/{param_set}/snp_density_heatmap.log",
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/plot_snp_density_heatmap.py"

rule coverage_density_heatmap:
    # All samples, one heatmap per (aligner, param_set): rows=sample,
    # columns=genomic window, color=mean depth.
    input:
        coverage_regions_for
    output:
        "results/mapping/{aligner}/{param_set}/metrics/coverage_density_heatmap.png",
    log:
        "logs/summary/{aligner}/{param_set}/coverage_density_heatmap.log",
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/plot_coverage_density_heatmap.py"
