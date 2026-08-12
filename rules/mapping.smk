# Align each sample once per (aligner, param_set) combination so mapping
# quality can be compared across parameter choices for the new prep.

rule bwa_mem_align:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
        idx=multiext(REF_FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam="results/mapping/bwa_mem/{param_set}/bam/{sample}.sorted.bam",
    params:
        extra=lambda wc: BWA_PARAMS[wc.param_set],
        ref=REF_FASTA,
    log:
        "logs/mapping/bwa_mem/{param_set}/{sample}.log",
    threads: get_mapping_threads()
    resources:
        mem_mb=get_mapping_mem_mb,
        runtime=get_mapping_runtime,
        partition=get_mapping_partition,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        bwa mem -t {threads} {params.extra} {params.ref} {input.r1} {input.r2} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam} -
        """

rule bowtie2_align:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
        idx=multiext(REF_FASTA, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    output:
        bam="results/mapping/bowtie2/{param_set}/bam/{sample}.sorted.bam",
    params:
        extra=lambda wc: BOWTIE2_PARAMS[wc.param_set],
        ref=REF_FASTA,
    log:
        "logs/mapping/bowtie2/{param_set}/{sample}.log",
    threads: get_mapping_threads()
    resources:
        mem_mb=get_mapping_mem_mb,
        runtime=get_mapping_runtime,
        partition=get_mapping_partition,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        bowtie2 -p {threads} {params.extra} -x {params.ref} -1 {input.r1} -2 {input.r2} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam} -
        """

rule index_bam:
    # CSI, not BAI: several wheat chromosomes (e.g. 3B ~830 Mbp) exceed
    # BAI's 2^29-1 (~512 Mbp) per-contig length limit.
    input:
        "results/mapping/{aligner}/{param_set}/bam/{sample}.sorted.bam",
    output:
        "results/mapping/{aligner}/{param_set}/bam/{sample}.sorted.bam.csi",
    conda:
        "../envs/qc.yml"
    shell:
        "samtools index -c {input}"

rule mark_duplicates:
    input:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.sorted.bam",
    output:
        bam="results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
        metrics="results/mapping/{aligner}/{param_set}/metrics/{sample}.dup_metrics.txt",
    log:
        "logs/mapping/{aligner}/{param_set}/{sample}.markdup.log",
    resources:
        mem_mb=get_mem_mb,
        runtime=get_runtime,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {output.metrics}) $(dirname {log})
        picard MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} &> {log}
        """

rule index_dedup_bam:
    # CSI, not BAI: several wheat chromosomes (e.g. 3B ~830 Mbp) exceed
    # BAI's 2^29-1 (~512 Mbp) per-contig length limit.
    input:
        "results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam",
    output:
        "results/mapping/{aligner}/{param_set}/bam/{sample}.dedup.bam.csi",
    conda:
        "../envs/qc.yml"
    shell:
        "samtools index -c {input}"
