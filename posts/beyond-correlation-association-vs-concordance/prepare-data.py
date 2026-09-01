"""Build the committed data file behind this post's figures and widget.

Provenance script. It is NOT executed when the site renders -- `index.qmd` contains
no code chunks and reads only the committed `methylation-subset.json`. Run it by
hand if the numbers in the post ever need regenerating:

    source .venv/bin/activate
    python posts/beyond-correlation-association-vs-concordance/prepare-data.py

Source: NCBI GEO accession GSE81211, "Genome-wide methylation approach identifies a
novel hypermethylated gene panel in ulcerative colitis". Illumina Infinium
HumanMethylation450 BeadChip (GPL13534), colonic mucosa: 3 normal, 8 ulcerative
colitis, plus the HCT116 colon cancer cell line. The series matrix is public and
de-identified; only the sample labels GEO itself publishes are carried through.

The 27 MB series matrix is downloaded to a temporary directory and never enters the
repository. Only the derived subset below is committed.
"""

import gzip
import json
import os
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ACCESSION = "GSE81211"
PLATFORM = "GPL13534"
URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/"
    f"{ACCESSION}/matrix/{ACCESSION}_series_matrix.txt.gz"
)

# Beta values on the 450K reach 0 and 1, where the M-value transform diverges.
# Clipping is the standard guard; the post says so explicitly.
CLIP = 0.001

N_SUBSET = 3000          # CpGs carried into the browser widget
SEED = 20260831          # fixed, so the committed subset is reproducible
HIST_BINS = 40

OUT = Path(__file__).parent / "methylation-subset.json"


def download(dest: Path) -> Path:
    """Fetch the series matrix into a scratch directory (not the repo)."""
    if dest.exists():
        print(f"using cached {dest}")
        return dest
    print(f"downloading {URL}")
    urllib.request.urlretrieve(URL, dest)
    return dest


def read_series_matrix(path: Path):
    """Return (sample titles keyed by GSM, beta DataFrame indexed by probe id)."""
    titles, gsms = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [f.strip('"') for f in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [f.strip('"') for f in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                break
        beta = pd.read_csv(fh, sep="\t", index_col=0, na_values=["null", "NA", ""])

    beta = beta[~beta.index.astype(str).str.startswith("!")]
    beta = beta.apply(pd.to_numeric, errors="coerce")
    return dict(zip(gsms, titles)), beta


def group_of(title: str) -> str:
    low = title.lower()
    if "normal" in low:
        return "normal"
    if "ulcerative" in low:
        return "uc"
    return "cancer-line"


def short_label(title: str) -> str:
    """'Ulcerative colitis 3' -> 'UC 3'; keeps axis labels readable."""
    return title.replace("Ulcerative colitis", "UC").replace(
        "Human colon cancer cell", "HCT116"
    )


def corr_matrices(mat: np.ndarray):
    """Pearson and Spearman across every sample pair. Columns are samples."""
    pearson = np.corrcoef(mat, rowvar=False)
    spearman = np.corrcoef(
        np.apply_along_axis(rankdata, 0, mat), rowvar=False
    )
    return pearson, spearman


def offdiag(m: np.ndarray) -> np.ndarray:
    return m[~np.eye(m.shape[0], dtype=bool)]


def summarise(name: str, m: np.ndarray, labels, groups) -> dict:
    """Print and return the spread of one coefficient matrix, split by group pair."""
    n = len(labels)
    mucosa = [i for i, g in enumerate(groups) if g in ("normal", "uc")]
    within_n, within_uc, across, vs_line = [], [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            v = m[i, j]
            if groups[i] == "cancer-line" or groups[j] == "cancer-line":
                vs_line.append(v)
            elif groups[i] == groups[j] == "normal":
                within_n.append(v)
            elif groups[i] == groups[j] == "uc":
                within_uc.append(v)
            else:
                across.append(v)

    mm = m[np.ix_(mucosa, mucosa)]
    out = {
        "min": float(offdiag(mm).min()),
        "max": float(offdiag(mm).max()),
        "mean": float(offdiag(mm).mean()),
        "normal_normal": float(np.mean(within_n)),
        "uc_uc": float(np.mean(within_uc)),
        "normal_uc": float(np.mean(across)),
        "vs_cancer_line": float(np.mean(vs_line)),
    }
    print(f"\n  {name} (11 mucosa samples, all pairs)")
    print(f"    range        {out['min']:.4f} - {out['max']:.4f}   mean {out['mean']:.4f}")
    print(f"    normal-normal {out['normal_normal']:.4f}")
    print(f"    UC-UC         {out['uc_uc']:.4f}")
    print(f"    normal-UC     {out['normal_uc']:.4f}")
    print(f"    vs HCT116     {out['vs_cancer_line']:.4f}")
    return out


def group_means(m: np.ndarray, groups) -> dict:
    """Mean coefficient for normal-normal, UC-UC and normal-UC sample pairs."""
    nn, uu, nu = [], [], []
    n = len(groups)
    for i in range(n):
        for j in range(i + 1, n):
            if "cancer-line" in (groups[i], groups[j]):
                continue
            bucket = nn if groups[i] == groups[j] == "normal" else (
                uu if groups[i] == groups[j] == "uc" else nu)
            bucket.append(m[i, j])
    same = (float(np.mean(nn)) + float(np.mean(uu))) / 2
    return {
        "normal_normal": float(np.mean(nn)),
        "uc_uc": float(np.mean(uu)),
        "normal_uc": float(np.mean(nu)),
        "separation": same - float(np.mean(nu)),
    }


def variability_sweep(beta: np.ndarray, sd: np.ndarray, groups, percents) -> list:
    """How the coefficients behave as the CpG set narrows to the most variable probes.

    This is the post's payoff: on all probes both coefficients are ~0.98 and neither
    separates disease groups. Filtering to variable CpGs pulls them apart, and the
    separation is wider under Spearman -- until the set gets so small that it inverts.
    """
    rows = []
    for pct in percents:
        keep = sd >= np.percentile(sd, 100 - pct)
        sub = beta[keep]
        p_mat, s_mat = corr_matrices(sub)
        rows.append({
            "percent": pct,
            "n_cpg": int(keep.sum()),
            "pearson": group_means(p_mat, groups),
            "spearman": group_means(s_mat, groups),
        })
        r = rows[-1]
        print(f"    top {pct:>3}%  n={r['n_cpg']:>6}  "
              f"Pearson sep {r['pearson']['separation']:+.4f}  "
              f"Spearman sep {r['spearman']['separation']:+.4f}")
    return rows


def main() -> None:
    scratch = Path(os.environ.get("TMPDIR", tempfile.gettempdir()))
    titles, beta_df = read_series_matrix(
        download(scratch / f"{ACCESSION}_series_matrix.txt.gz")
    )

    n_raw = len(beta_df)
    beta_df = beta_df.dropna()
    print(f"{ACCESSION}: {n_raw} probes, {len(beta_df)} complete across "
          f"{beta_df.shape[1]} samples")

    labels = [short_label(titles[g]) for g in beta_df.columns]
    groups = [group_of(titles[g]) for g in beta_df.columns]
    for gsm, lab, grp in zip(beta_df.columns, labels, groups):
        print(f"  {gsm}  {lab:<10} {grp}")

    beta = beta_df.to_numpy(dtype=float)
    clipped = np.clip(beta, CLIP, 1 - CLIP)
    mval = np.log2(clipped / (1 - clipped))

    print("\n=== full data, all %d CpGs ===" % len(beta_df))
    b_pearson, b_spearman = corr_matrices(beta)
    m_pearson, m_spearman = corr_matrices(mval)

    stats = {
        "beta_pearson": summarise("beta / Pearson", b_pearson, labels, groups),
        "beta_spearman": summarise("beta / Spearman", b_spearman, labels, groups),
        "mval_pearson": summarise("M-value / Pearson", m_pearson, labels, groups),
        "mval_spearman": summarise("M-value / Spearman", m_spearman, labels, groups),
    }

    # The post's central claim: a monotone transform leaves Spearman untouched.
    # The Pearson shift is reported two ways, because pooling them is misleading:
    # the cancer-line pairs move an order of magnitude more than the mucosa pairs,
    # so an all-pairs mean is dominated by the eleven comparisons against HCT116.
    drift = np.abs(b_spearman - m_spearman).max()
    delta = np.abs(b_pearson - m_pearson)
    n_s = len(groups)
    muc_shift, line_shift = [], []
    for i in range(n_s):
        for j in range(i + 1, n_s):
            target = line_shift if "cancer-line" in (groups[i], groups[j]) else muc_shift
            target.append(delta[i, j])
    shift_mucosa = float(np.mean(muc_shift))
    shift_line = float(np.mean(line_shift))
    shift_all = float(np.abs(offdiag(b_pearson) - offdiag(m_pearson)).mean())
    print(f"\n  Spearman drift beta -> M-value: max {drift:.2e}")
    print(f"  Pearson shift, mucosa pairs:    {shift_mucosa:.4f} "
          f"(max {max(muc_shift):.4f})")
    print(f"  Pearson shift, cancer-line:     {shift_line:.4f}")
    print(f"  Pearson shift, all pairs:       {shift_all:.4f}")

    # Shape of the beta distribution, pooled across mucosa samples.
    mucosa_cols = [i for i, g in enumerate(groups) if g != "cancer-line"]
    counts, edges = np.histogram(
        beta[:, mucosa_cols].ravel(), bins=HIST_BINS, range=(0, 1)
    )
    frac = counts / counts.sum()
    ends = float(frac[:HIST_BINS // 10].sum() + frac[-HIST_BINS // 10:].sum())
    middle = float(frac[HIST_BINS * 3 // 10: HIST_BINS * 7 // 10].sum())
    print(f"\n  beta below 0.1 or above 0.9: {ends:.3f}")
    print(f"  beta between 0.3 and 0.7:    {middle:.3f}")

    # Widget subset: a fixed-seed random sample of probes, with each probe's
    # spread across the mucosa samples so the variability filter is meaningful.
    sd_all = beta[:, mucosa_cols].std(axis=1, ddof=1)
    mucosa_groups = [groups[i] for i in mucosa_cols]
    print("\n  narrowing to the most variable CpGs")
    sweep = variability_sweep(
        beta[:, mucosa_cols], sd_all, mucosa_groups, [100, 50, 25, 10, 5, 2, 1]
    )

    rng = np.random.default_rng(SEED)
    idx = np.sort(rng.choice(len(beta_df), size=N_SUBSET, replace=False))
    sd = sd_all[idx]

    payload = {
        "accession": ACCESSION,
        "platform": PLATFORM,
        "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={ACCESSION}",
        "n_cpg_total": int(len(beta_df)),
        "n_cpg_subset": int(N_SUBSET),
        "clip": CLIP,
        "seed": SEED,
        "samples": [
            {"gsm": g, "label": lab, "group": grp}
            for g, lab, grp in zip(beta_df.columns, labels, groups)
        ],
        "cpgs": [str(c) for c in beta_df.index[idx]],
        # Encoded as integers to keep the committed file small: beta x1000 and
        # sd x10000. The browser divides them back out on load.
        "beta_scale": 1000,
        "sd_scale": 10000,
        "beta": np.rint(beta[idx] * 1000).astype(int).tolist(),
        "sd": np.rint(sd * 10000).astype(int).tolist(),
        "histogram": {
            "edges": np.round(edges, 3).tolist(),
            "fraction": np.round(frac, 5).tolist(),
            "frac_at_ends": round(ends, 4),
            "frac_in_middle": round(middle, 4),
        },
        "full_stats": stats,
        "variability_sweep": sweep,
        "spearman_drift_beta_to_mval": float(drift),
        "pearson_shift_beta_to_mval_mucosa": shift_mucosa,
        "pearson_shift_beta_to_mval_cancer_line": shift_line,
        "pearson_shift_beta_to_mval_all_pairs": shift_all,
        "matrices": {
            "labels": labels,
            "groups": groups,
            "beta_pearson": np.round(b_pearson, 4).tolist(),
            "beta_spearman": np.round(b_spearman, 4).tolist(),
            "mval_pearson": np.round(m_pearson, 4).tolist(),
            "mval_spearman": np.round(m_spearman, 4).tolist(),
        },
    }

    OUT.write_text(json.dumps(payload, separators=(",", ":")))
    print(f"\nwrote {OUT} ({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
