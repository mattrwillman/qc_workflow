import numpy as np
import pandas as pd

# mosdepth's *.mosdepth.global.dist.txt gives, per chrom (+ "total"), the
# proportion of bases with coverage >= depth for every observed depth —
# i.e. a survival function, not a raw depth list. Convert it to a
# probability mass function to get mean/stdev/percentiles without ever
# touching a per-base file (skipped via --no-per-base).
df = pd.read_csv(
    snakemake.input["dist"], sep="\t", names=["chrom", "depth", "proportion"]
)
df = df[df["chrom"] == "total"].sort_values("depth").reset_index(drop=True)

depths = df["depth"].to_numpy(dtype=float)
survival = df["proportion"].to_numpy(dtype=float)

# pmf[i] = P(X == depths[i]) = survival(depths[i]) - survival(depths[i]+1);
# the last (max-depth) row has no "+1" row, so its remaining mass is just
# its own survival value.
pmf = -np.diff(survival)
pmf = np.append(pmf, survival[-1])
pmf = np.clip(pmf, 0, None)
pmf = pmf / pmf.sum()

mean = float(np.sum(depths * pmf))
variance = float(np.sum(((depths - mean) ** 2) * pmf))
stdev = variance**0.5
cv = stdev / mean if mean > 0 else float("nan")

cdf = np.cumsum(pmf)


def percentile(p):
    idx = int(np.searchsorted(cdf, p))
    idx = min(idx, len(depths) - 1)
    return depths[idx]


row = {
    "sample": snakemake.wildcards.sample,
    "aligner": snakemake.wildcards.aligner,
    "param_set": snakemake.wildcards.param_set,
    "mean": round(mean, 4),
    "stdev": round(stdev, 4),
    "cv": round(cv, 4),
    "p10": percentile(0.10),
    "p25": percentile(0.25),
    "median": percentile(0.50),
    "p75": percentile(0.75),
    "p90": percentile(0.90),
}

pd.DataFrame([row]).to_csv(snakemake.output[0], sep="\t", index=False)
