# qc_workflow

FastQC + mapping-parameter QC sweep for evaluating a new tagmentation
library prep.

## What it does

1. **FastQC** on raw R1/R2 reads for every sample.
2. **Alignment** of every sample against the configured reference, once per
   `(aligner, param_set)` combination defined in `config/config.yml`
   (`bwa_mem_params` / `bowtie2_params`), so mapping quality can be compared
   across parameter choices.
3. Per alignment: **duplicate marking** (Picard), **flagstat/samtools stats**,
   **insert size distribution** (Picard), and **coverage/depth** (mosdepth) —
   the four signals that matter most for judging a tagmentation prep
   (fragment size distribution, duplication/over-amplification rate,
   mapping rate, and coverage uniformity).
4. **SNP calling**: `bcftools mpileup | bcftools call` + `bcftools stats` per
   `(aligner, param_set, sample)`, plus:
   - a per-sample **SNP density plot** (SNPs per window, one track per
     chromosome/contig), and
   - a **cross-sample SNP density heatmap**, one per `(aligner, param_set)`
     — rows are samples, columns are genomic windows, color is SNP count,
     so uneven regions/outlier samples are visible at a glance.
5. **Read coverage density**: windowed mean depth from mosdepth (`--by`,
   same window size as the SNP density plots), with the same per-sample
   line plot + cross-sample heatmap pairing as the SNP density above.
6. **Cross-sample coverage summary** (`results/summary/coverage_summary.tsv`):
   mean, stdev, CV, and percentiles (p10/p25/median/p75/p90), one row per
   `(sample, aligner, param_set)`. Derived from mosdepth's cumulative depth
   distribution rather than min/mean/max — for skim-seq coverage well below
   1x, min is trivially always 0 and breadth-at-depth thresholds (e.g. %
   >=1x) aren't informative either; CV captures uniformity and percentiles
   capture the actual shape of a mostly-sub-1x distribution.
7. **MultiQC** aggregates everything (FastQC, mapping metrics, bcftools
   stats) into `results/multiqc/qc_report.html`.

The SNP density and coverage density windows are both controlled by
`density_window_bp` in `config/config.yml` (default 1 Mb) so the two line
up on the same genomic bins.

> **Caveat on the SNP data**: with a single library per genotype and no
> replicates, SNP counts/het ratios are confounded by each genotype's
> evolutionary distance to the reference — a "noisy" sample could mean a
> divergent genotype rather than a library-prep problem. Treat SNP density
> as a secondary/exploratory signal; FastQC and the alignment-based metrics
> above (duplication rate, insert size, coverage uniformity, mapping rate)
> are the more reliable read on the tagmentation prep itself.

## Setup

Two conda envs are involved:

- `environment.yml` — the **controller** env (Snakemake itself + the
  `cluster-generic` executor plugin needed for `profiles/slurm`). Activate
  this before invoking `snakemake`.
- `envs/qc.yml` — **per-rule tool** env (fastqc, bwa, bowtie2, samtools,
  bcftools, picard, mosdepth, multiqc, numpy/pandas/matplotlib for the
  density plots), pulled in automatically via `--use-conda`. You don't
  activate this one yourself.

```bash
conda env create -f environment.yml
conda activate qc_workflow
```

Then edit `config/config.yml`:

- `fastq_dir`, `r1_suffix`, `r2_suffix` — where raw fastqs live and how R1/R2
  are named. Samples are auto-discovered from this directory.
- `reference.fasta` — path to the reference genome to map against.
- `aligners`, `bwa_mem_params`, `bowtie2_params` — which parameter sets to
  sweep. Add/remove entries to test more or fewer combinations.

## Run

From inside `qc_workflow/`, with the `qc_workflow` controller env active:

```bash
snakemake -s Snakefile --configfile config/config.yml \
    --use-conda --profile profiles/local
```

Swap `profiles/local` for `profiles/slurm` to run on Atlas. The controller
must run on the login node inside `tmux`, not via `sbatch` — the SLURM
executor plugin fails silently when run from a compute node because it
lacks sacct/squeue access needed for job monitoring.

## Output layout

```
results/
  fastqc/raw/{sample}_{R1,R2}_fastqc.{html,zip}
  mapping/{aligner}/{param_set}/bam/{sample}.{sorted,dedup}.bam
  mapping/{aligner}/{param_set}/metrics/{sample}.{flagstat,samtools_stats,
      insert_size_metrics,dup_metrics,mosdepth.summary,
      mosdepth.global.dist,coverage_stats}.txt
  mapping/{aligner}/{param_set}/metrics/{sample}.regions.bed.gz
  mapping/{aligner}/{param_set}/metrics/{sample}.coverage_density.png
  mapping/{aligner}/{param_set}/metrics/coverage_density_heatmap.png
  mapping/{aligner}/{param_set}/variants/{sample}.vcf.gz
  mapping/{aligner}/{param_set}/variants/{sample}.bcftools_stats.txt
  mapping/{aligner}/{param_set}/variants/{sample}.snp_density.{tsv,png}
  mapping/{aligner}/{param_set}/variants/snp_density_heatmap.png
  summary/coverage_summary.tsv
  multiqc/qc_report.html
```

## Development notes

Developed with [Claude Code](https://claude.com/claude-code) (Anthropic),
under the workflow author's direction and review.
