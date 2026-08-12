import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Cap the number of chromosomes/contigs plotted so a scaffold-heavy
# assembly doesn't blow up into hundreds of subplots.
MAX_CHROMS = 30

df = pd.read_csv(
    snakemake.input[0], sep="\t", names=["chrom", "window_start", "snp_count"]
)

top_chroms = (
    df.groupby("chrom")["snp_count"].sum().sort_values(ascending=False).head(MAX_CHROMS)
)
chroms = list(top_chroms.index)
df = df[df["chrom"].isin(chroms)]

fig, axes = plt.subplots(len(chroms), 1, figsize=(10, 1.4 * len(chroms)), sharex=False)
if len(chroms) == 1:
    axes = [axes]

for ax, chrom in zip(axes, chroms):
    sub = df[df["chrom"] == chrom].sort_values("window_start")
    ax.plot(sub["window_start"] / 1e6, sub["snp_count"], lw=0.8)
    ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=8)

axes[-1].set_xlabel("Position (Mb)")
fig.suptitle(
    f"SNP density — {snakemake.wildcards.sample} "
    f"({snakemake.wildcards.aligner}/{snakemake.wildcards.param_set})"
)
fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=150)
