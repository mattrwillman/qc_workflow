import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Heatmap rows are one-per-sample, so keep fewer chromosome tracks than the
# per-sample line plot to keep the figure a manageable size.
MAX_CHROMS = 15


def sample_name(path):
    return re.sub(r"\.snp_density\.tsv$", "", path.split("/")[-1])


frames = []
for path in snakemake.input:
    df = pd.read_csv(path, sep="\t", names=["chrom", "window_start", "snp_count"])
    df["sample"] = sample_name(path)
    frames.append(df)
all_df = pd.concat(frames, ignore_index=True)

top_chroms = (
    all_df.groupby("chrom")["snp_count"].sum().sort_values(ascending=False).head(MAX_CHROMS)
)
chroms = list(top_chroms.index)
all_df = all_df[all_df["chrom"].isin(chroms)]
samples = sorted(all_df["sample"].unique())

fig, axes = plt.subplots(
    len(chroms), 1, figsize=(12, 0.35 * len(samples) * len(chroms) + len(chroms)), sharex=False
)
if len(chroms) == 1:
    axes = [axes]

# Robust shared color scale across all chromosomes/samples — a single
# outlier window shouldn't wash out the rest of the heatmap.
vmax = max(float(all_df["snp_count"].quantile(0.99)), 1.0)

im = None
for ax, chrom in zip(axes, chroms):
    sub = all_df[all_df["chrom"] == chrom]
    pivot = sub.pivot_table(
        index="sample", columns="window_start", values="snp_count", fill_value=0
    )
    pivot = pivot.reindex(index=samples).reindex(columns=sorted(pivot.columns))
    im = ax.imshow(
        pivot.values,
        aspect="auto",
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        extent=[pivot.columns.min() / 1e6, pivot.columns.max() / 1e6, len(samples), 0],
    )
    ax.set_yticks(np.arange(len(samples)) + 0.5)
    ax.set_yticklabels(samples, fontsize=7)
    ax.set_ylabel(chrom, fontsize=9)

axes[-1].set_xlabel("Position (Mb)")
fig.suptitle(
    f"SNP density across samples ({snakemake.wildcards.aligner}/{snakemake.wildcards.param_set})"
)
fig.colorbar(im, ax=axes, shrink=0.5, label="SNPs per window")
fig.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
