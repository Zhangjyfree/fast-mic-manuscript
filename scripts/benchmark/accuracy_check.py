#!/usr/bin/env python3
"""
Correctness check: compare fast-mic vs COBRApy growth rates.
Outputs a stats TSV and a scatter-ready TSV for plotting.

Usage:
    python accuracy_check.py <fastmic.tsv> <cobra.tsv> <stats_out.tsv> <scatter_out.tsv>
"""

import sys
import math


def parse_tsv(path):
    rows = {}
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("model_id\t"):
                header = line.split("\t")
                continue
            if header is None:
                continue
            fields = line.split("\t")
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            mid = row["model_id"]
            rows[mid] = float(row.get("growth_rate", 0) or 0)
    return rows


def pearson_r(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    dx = math.sqrt(sum((xi - mx) ** 2 for xi in x))
    dy = math.sqrt(sum((yi - my) ** 2 for yi in y))
    return num / (dx * dy) if dx * dy > 0 else float("nan")


def main():
    if len(sys.argv) != 5:
        print(
            "Usage: accuracy_check.py <fastmic.tsv> <cobra.tsv> "
            "<stats_out.tsv> <scatter_out.tsv>",
            file=sys.stderr,
        )
        sys.exit(1)

    fastmic = parse_tsv(sys.argv[1])
    cobra   = parse_tsv(sys.argv[2])
    stats_path   = sys.argv[3]
    scatter_path = sys.argv[4]

    common = sorted(set(fastmic) & set(cobra))
    if not common:
        print("ERROR: no models in common between the two TSVs", file=sys.stderr)
        sys.exit(1)

    fm = [fastmic[m] for m in common]
    cb = [cobra[m]   for m in common]

    # ── Statistics ──
    n = len(common)
    mae  = sum(abs(f - c) for f, c in zip(fm, cb)) / n
    rmse = math.sqrt(sum((f - c) ** 2 for f, c in zip(fm, cb)) / n)
    r    = pearson_r(fm, cb)
    r2   = r ** 2 if not math.isnan(r) else float("nan")

    # % with growth > 0 in both (non-zero pairs)
    nonzero = [(f, c) for f, c in zip(fm, cb) if f > 1e-6 and c > 1e-6]
    pct_agree_1 = sum(
        1 for f, c in nonzero if abs(f - c) / c < 0.01
    ) / max(len(nonzero), 1) * 100
    pct_agree_5 = sum(
        1 for f, c in nonzero if abs(f - c) / c < 0.05
    ) / max(len(nonzero), 1) * 100

    n_both_zero  = sum(1 for f, c in zip(fm, cb) if f < 1e-6 and c < 1e-6)
    n_fm_only    = sum(1 for f, c in zip(fm, cb) if f > 1e-6 and c < 1e-6)
    n_cobra_only = sum(1 for f, c in zip(fm, cb) if f < 1e-6 and c > 1e-6)
    n_both_pos   = len(nonzero)

    # Print summary to stderr
    print(f"\n{'='*52}", file=sys.stderr)
    print(f"  Correctness: fast-mic vs COBRApy ({n} models)", file=sys.stderr)
    print(f"{'='*52}", file=sys.stderr)
    print(f"  Pearson r:          {r:.6f}", file=sys.stderr)
    print(f"  R²:                 {r2:.6f}", file=sys.stderr)
    print(f"  MAE (h⁻¹):         {mae:.2e}", file=sys.stderr)
    print(f"  RMSE (h⁻¹):        {rmse:.2e}", file=sys.stderr)
    print(f"  Both grow (>0):     {n_both_pos}/{n}", file=sys.stderr)
    print(f"  Agree within 1%:    {pct_agree_1:.1f}% of growing pairs", file=sys.stderr)
    print(f"  Agree within 5%:    {pct_agree_5:.1f}% of growing pairs", file=sys.stderr)
    print(f"  Both zero:          {n_both_zero}", file=sys.stderr)
    print(f"  fast-mic only > 0:  {n_fm_only}", file=sys.stderr)
    print(f"  COBRApy only > 0:   {n_cobra_only}", file=sys.stderr)
    print(f"{'='*52}\n", file=sys.stderr)

    # ── Write stats TSV ──
    with open(stats_path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"n_models\t{n}\n")
        fh.write(f"pearson_r\t{r:.8f}\n")
        fh.write(f"r_squared\t{r2:.8f}\n")
        fh.write(f"mae\t{mae:.8e}\n")
        fh.write(f"rmse\t{rmse:.8e}\n")
        fh.write(f"n_both_positive\t{n_both_pos}\n")
        fh.write(f"pct_agree_1pct\t{pct_agree_1:.4f}\n")
        fh.write(f"pct_agree_5pct\t{pct_agree_5:.4f}\n")
        fh.write(f"n_both_zero\t{n_both_zero}\n")
        fh.write(f"n_fastmic_only_positive\t{n_fm_only}\n")
        fh.write(f"n_cobra_only_positive\t{n_cobra_only}\n")

    # ── Write scatter TSV ──
    with open(scatter_path, "w") as fh:
        fh.write("model_id\tfastmic_growth\tcobra_growth\tabs_diff\trel_diff_pct\n")
        for mid, f, c in zip(common, fm, cb):
            abs_diff = abs(f - c)
            rel_diff = abs(f - c) / c * 100 if c > 1e-6 else float("nan")
            rel_str  = f"{rel_diff:.4f}" if not math.isnan(rel_diff) else "NA"
            fh.write(f"{mid}\t{f:.8f}\t{c:.8f}\t{abs_diff:.2e}\t{rel_str}\n")

    print(f"Stats:   {stats_path}", file=sys.stderr)
    print(f"Scatter: {scatter_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
