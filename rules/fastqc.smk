# Raw-read quality: FastQC per sample per mate.

rule fastqc_raw:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
    output:
        html_r1="results/fastqc/raw/{sample}_R1_fastqc.html",
        html_r2="results/fastqc/raw/{sample}_R2_fastqc.html",
        zip_r1="results/fastqc/raw/{sample}_R1_fastqc.zip",
        zip_r2="results/fastqc/raw/{sample}_R2_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    threads: 2
    resources:
        mem_mb=get_mem_mb,
        runtime=get_runtime,
    conda:
        "../envs/qc.yml"
    shell:
        """
        mkdir -p results/fastqc/raw logs/fastqc
        tmpdir=$(mktemp -d)
        fastqc --threads {threads} --outdir "$tmpdir" {input.r1} {input.r2} &> {log}
        mv "$tmpdir"/*_R1*_fastqc.html {output.html_r1}
        mv "$tmpdir"/*_R1*_fastqc.zip {output.zip_r1}
        mv "$tmpdir"/*_R2*_fastqc.html {output.html_r2}
        mv "$tmpdir"/*_R2*_fastqc.zip {output.zip_r2}
        rm -rf "$tmpdir"
        """
