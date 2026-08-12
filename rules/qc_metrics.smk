# Mapping-quality metrics computed on the deduplicated BAM for each
# (aligner, param_set, sample) combination.

rule samtools_flagstat:
    input:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
        bai="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam.csi",
    output:
        "results/mapping/{aligner}/{param_set}/metrics/{sample}.flagstat.txt",
    conda:
        "../envs/qc.yml"
    shell:
        "samtools flagstat {input.bam} > {output}"

rule samtools_stats:
    input:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
        bai="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam.csi",
    output:
        "results/mapping/{aligner}/{param_set}/metrics/{sample}.samtools_stats.txt",
    conda:
        "../envs/qc.yml"
    shell:
        "samtools stats {input.bam} > {output}"

rule picard_insert_size:
    input:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
        bai="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam.csi",
    output:
        metrics="results/mapping/{aligner}/{param_set}/metrics/{sample}.insert_size_metrics.txt",
        hist="results/mapping/{aligner}/{param_set}/metrics/{sample}.insert_size_histogram.pdf",
    log:
        "logs/qc_metrics/{aligner}/{param_set}/{sample}.insert_size.log",
    resources:
        mem_mb=get_mem_mb,
        runtime=get_runtime,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p $(dirname {output.metrics}) $(dirname {log})
        picard CollectInsertSizeMetrics \
            INPUT={input.bam} \
            OUTPUT={output.metrics} \
            HISTOGRAM_FILE={output.hist} &> {log}
        """

rule mosdepth_coverage:
    input:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
        bai="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam.csi",
    output:
        summary="results/mapping/{aligner}/{param_set}/metrics/{sample}.mosdepth.summary.txt",
        dist="results/mapping/{aligner}/{param_set}/metrics/{sample}.mosdepth.global.dist.txt",
        # Windowed mean depth (--by, no per-base file) for the coverage
        # density plots — same window size as the SNP density plots.
        regions="results/mapping/{aligner}/{param_set}/metrics/{sample}.regions.bed.gz",
    params:
        prefix="results/mapping/{aligner}/{param_set}/metrics/{sample}",
        window=DENSITY_WINDOW,
    threads: 4
    resources:
        mem_mb=get_mem_mb,
        runtime=get_runtime,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p $(dirname {output.summary})
        mosdepth --threads {threads} --no-per-base --by {params.window} {params.prefix} {input.bam}
        """

rule coverage_density_plot:
    # Windowed mean-depth line plot, one track per chromosome — analogous
    # to the SNP density plot but for read coverage.
    input:
        "results/mapping/{aligner}/{param_set}/metrics/{sample}.regions.bed.gz",
    output:
        "results/mapping/{aligner}/{param_set}/metrics/{sample}.coverage_density.png",
    log:
        "logs/qc_metrics/{aligner}/{param_set}/{sample}.coverage_density_plot.log",
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/plot_coverage_density.py"

rule coverage_stats:
    # Mean/stdev/CV + percentiles derived from mosdepth's cumulative depth
    # distribution. Skim-seq coverage is well below 1x, so breadth-at-depth
    # thresholds (e.g. % >=1x) aren't informative — percentiles and CV are.
    input:
        dist="results/mapping/{aligner}/{param_set}/metrics/{sample}.mosdepth.global.dist.txt",
    output:
        "results/mapping/{aligner}/{param_set}/metrics/{sample}.coverage_stats.tsv",
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/coverage_stats.py"


# ── SNP calling ──────────────────────────────────────────────────────────────
# Single-library-per-genotype SNP counts are confounded by each genotype's
# relatedness to the reference, so treat these as a secondary, exploratory
# signal alongside (not instead of) the alignment-based QC metrics above.

rule bcftools_call:
    input:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
        bai="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam.csi",
        fai=REF_FASTA + ".fai",
    output:
        vcf="results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz",
    params:
        ref=REF_FASTA,
    log:
        "logs/qc_metrics/{aligner}/{param_set}/{sample}.bcftools_call.log",
    threads: 4
    resources:
        mem_mb=get_mem_mb,
        runtime=get_runtime,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        bcftools mpileup --threads {threads} -f {params.ref} {input.bam} 2> {log} \
            | bcftools call --threads {threads} -mv -Oz -o {output.vcf} 2>> {log}
        """

rule bcftools_index:
    # CSI, not TBI: several wheat chromosomes (e.g. 3B ~830 Mbp) exceed
    # TBI's 2^29-1 (~512 Mbp) per-contig length limit, same issue as the
    # BAM indexing rules in mapping.smk.
    input:
        "results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz",
    output:
        "results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz.csi",
    conda:
        "../envs/qc.yml"
    shell:
        "bcftools index -c {input}"

rule bcftools_stats:
    input:
        vcf="results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz",
        csi="results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz.csi",
    output:
        "results/mapping/{aligner}/{param_set}/variants/{sample}.bcftools_stats.txt",
    conda:
        "../envs/qc.yml"
    shell:
        "bcftools stats {input.vcf} > {output}"

rule snp_density_windows:
    input:
        vcf="results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz",
        csi="results/mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz.csi",
    output:
        "results/mapping/{aligner}/{param_set}/variants/{sample}.snp_density.tsv",
    params:
        window=DENSITY_WINDOW,
    conda:
        "../envs/qc.yml"
    shell:
        r"""
        bcftools query -f '%CHROM\t%POS\n' {input.vcf} \
            | awk -v w={params.window} 'BEGIN{{OFS="\t"}}{{bin=int(($2-1)/w)*w; print $1,bin}}' \
            | sort -k1,1 -k2,2n \
            | uniq -c \
            | awk 'BEGIN{{OFS="\t"}}{{print $2,$3,$1}}' > {output}
        """

rule snp_density_plot:
    input:
        "results/mapping/{aligner}/{param_set}/variants/{sample}.snp_density.tsv",
    output:
        "results/mapping/{aligner}/{param_set}/variants/{sample}.snp_density.png",
    log:
        "logs/qc_metrics/{aligner}/{param_set}/{sample}.snp_density_plot.log",
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/plot_snp_density.py"
