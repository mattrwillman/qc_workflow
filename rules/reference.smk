# Reference genome indexing, shared across all aligners/param sets.

rule bwa_index:
    input:
        REF_FASTA,
    output:
        multiext(REF_FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/reference/bwa_index.log",
    resources:
        mem_mb=get_index_mem_mb,
        runtime=get_index_runtime,
        partition=get_index_partition,
    conda:
        "../envs/qc.yml"
    shell:
        "bwa index {input} &> {log}"

rule bowtie2_index:
    input:
        REF_FASTA,
    output:
        multiext(REF_FASTA, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log:
        "logs/reference/bowtie2_index.log",
    threads: get_index_threads()
    resources:
        mem_mb=get_index_mem_mb,
        runtime=get_index_runtime,
        partition=get_index_partition,
    conda:
        "../envs/qc.yml"
    shell:
        "bowtie2-build --threads {threads} --large-index {input} {input} &> {log}"

rule samtools_faidx:
    input:
        REF_FASTA,
    output:
        REF_FASTA + ".fai",
    log:
        "logs/reference/samtools_faidx.log",
    conda:
        "../envs/qc.yml"
    shell:
        "samtools faidx {input} &> {log}"
