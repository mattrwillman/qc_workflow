# Aggregate FastQC + mapping metrics across all samples, aligners, and
# parameter sets into a single browsable report.

rule multiqc:
    input:
        expand(
            "results/fastqc/raw/{sample}_{read}_fastqc.zip",
            sample=SAMPLES, read=["R1", "R2"]
        ),
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.flagstat.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.samtools_stats.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.insert_size_metrics.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.dup_metrics.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.mosdepth.summary.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/variants/{sample}.bcftools_stats.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
    output:
        "results/multiqc/qc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/qc.yml"
    shell:
        """
        multiqc results/fastqc results/mapping \
            --filename qc_report.html \
            --outdir results/multiqc \
            --force &> {log}
        """
