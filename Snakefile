# qc_workflow — FastQC + mapping-parameter QC sweep for a new tagmentation
# library prep.
#
# Run with (from inside this directory):
#   snakemake -s Snakefile --configfile config/config.yml \
#       --use-conda --profile profiles/local
#
# Samples are discovered by scanning `fastq_dir` for files matching
# r1_suffix/r2_suffix. Each sample is mapped once per (aligner, param_set)
# combination in bwa_mem_params / bowtie2_params, so mapping-quality metrics
# can be compared across parameter choices for the new prep.

import os

configfile: "config/config.yml"

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# ── globals ──────────────────────────────────────────────────────────────────
FASTQ_DIR  = config["fastq_dir"]
R1_SUFFIX  = config["r1_suffix"]
R2_SUFFIX  = config["r2_suffix"]
REF_FASTA  = config["reference"]["fasta"]

ALIGNERS        = config["aligners"]
BWA_PARAMS      = config.get("bwa_mem_params", {})
BOWTIE2_PARAMS  = config.get("bowtie2_params", {})
DENSITY_WINDOW  = int(config.get("density_window_bp", 1000000))

SAMPLES, = glob_wildcards(os.path.join(FASTQ_DIR, "{sample}" + R1_SUFFIX))
if not SAMPLES:
    raise WorkflowError(
        f"No samples found in {FASTQ_DIR} matching *{R1_SUFFIX} / *{R2_SUFFIX}"
    )

def aligner_param_sets(aligner):
    return BWA_PARAMS if aligner == "bwa_mem" else BOWTIE2_PARAMS

# (aligner, param_set) pairs actually requested in config
ALIGNER_PARAM_PAIRS = [
    (aligner, param_set)
    for aligner in ALIGNERS
    for param_set in aligner_param_sets(aligner)
]


# ── helper functions ────────────────────────────────────────────────────────
def fastq_r1(wildcards):
    return os.path.join(FASTQ_DIR, wildcards.sample + R1_SUFFIX)

def fastq_r2(wildcards):
    return os.path.join(FASTQ_DIR, wildcards.sample + R2_SUFFIX)

def mapping_args(wildcards):
    return aligner_param_sets(wildcards.aligner)[wildcards.param_set]

def get_profile(name="default"):
    return config["resource_profiles"].get(name, config["resource_profiles"]["default"])

def get_threads(wc=None):
    return int(get_profile()["threads"])

def get_mem_mb(wc=None):
    return int(get_profile()["mem_mb"])

def get_runtime(wc=None):
    return int(get_profile()["runtime"])

def get_partition(wc=None):
    return get_profile()["partition"]

def get_index_threads(wc=None):
    return int(get_profile("genome_index")["threads"])

def get_index_mem_mb(wc=None):
    return int(get_profile("genome_index")["mem_mb"])

def get_index_runtime(wc=None):
    return int(get_profile("genome_index")["runtime"])

def get_index_partition(wc=None):
    return get_profile("genome_index")["partition"]

def get_mapping_threads(wc=None):
    return int(get_profile("mapping")["threads"])

def get_mapping_mem_mb(wc=None):
    return int(get_profile("mapping")["mem_mb"])

def get_mapping_runtime(wc=None):
    return int(get_profile("mapping")["runtime"])

def get_mapping_partition(wc=None):
    return get_profile("mapping")["partition"]

def snp_density_tsvs_for(wc):
    return [
        f"results/mapping/{wc.aligner}/{wc.param_set}/variants/{sample}.snp_density.tsv"
        for sample in SAMPLES
    ]

def coverage_regions_for(wc):
    return [
        f"results/mapping/{wc.aligner}/{wc.param_set}/metrics/{sample}.regions.bed.gz"
        for sample in SAMPLES
    ]


# ── lifecycle ────────────────────────────────────────────────────────────────
onstart:
    os.makedirs("logs/slurm", exist_ok=True)


# ── default target ──────────────────────────────────────────────────────────
rule all:
    input:
        expand(
            "results/fastqc/raw/{sample}_{read}_fastqc.html",
            sample=SAMPLES, read=["R1", "R2"]
        ),
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.flagstat.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.insert_size_metrics.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.mosdepth.summary.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.dup_metrics.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/variants/{sample}.bcftools_stats.txt"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/variants/{sample}.snp_density.png"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/{sample}.coverage_density.png"
            for sample in SAMPLES
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/variants/snp_density_heatmap.png"
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        [
            f"results/mapping/{aligner}/{param_set}/metrics/coverage_density_heatmap.png"
            for aligner, param_set in ALIGNER_PARAM_PAIRS
        ],
        "results/summary/coverage_summary.tsv",
        "results/multiqc/qc_report.html"


# ── rule modules ─────────────────────────────────────────────────────────────
include: "rules/fastqc.smk"
include: "rules/reference.smk"
include: "rules/mapping.smk"
include: "rules/qc_metrics.smk"
include: "rules/summary.smk"
include: "rules/multiqc.smk"
