#!/usr/bin/env python
"""Strain-level statistics for the generalist-vs-specialist contrast (Fig S4, Table S9).

Answers Reviewer 1, Comment 6 ("statistical inference largely absent") by treating the
probiotic STRAIN, not the pair, as the unit of replication: Akkermansia n = 6,
Lactobacillus n = 10. Every number quoted in the Results and in the Fig S4A annotation is
computed here, so nothing in the figure is hard-coded.

What it computes
  1. Per-strain mutualism % at L5 (pectin peak) and L6 (free glucose), over viable pairs.
  2. Genus mean with a 95% bootstrap CI (resampling STRAINS within a genus).
  3. Generalist vs specialist at L5:
       - ratio of genus means and its 95% bootstrap CI
       - difference of genus means (percentage points) and its 95% bootstrap CI
       - Mann-Whitney U (two-sided, exact) across strains
  4. The paired L5 -> L6 decline within each genus: Wilcoxon signed-rank (exact) and a
     paired bootstrap CI of the mean drop.

Bootstrap detail: strains are resampled with replacement WITHIN each genus, B = 10,000,
percentile CIs, fixed seed (0) so the numbers are reproducible. A viable pair = both
members grow in monoculture (> 1e-4 / h); mutualism % is taken over viable pairs.

Outputs (results/figS4/)
  per_strain_mutualism_L5_L6.tsv   system, strain, L5_mut_pct, L6_mut_pct
  genus_ci_summary.tsv             system, level, mean, lo, hi, n
  strain_level_tests.tsv           quantity, estimate, ci_lo, ci_hi, test, statistic,
                                   p_value, n_akk, n_lac, note

Usage: python3 scripts/stats_strain_level.py
"""
import csv, os
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon

FM = "."
VIABLE = 1e-4
B = 10_000
SEED = 0
LEVELS = {"L5": "L5_pectin", "L6": "L6_resistant_starch"}
SYS = {"akk": ("Akkermansia", f"{FM}/results/akk_vs_uhgg"),
       "lac": ("Lactobacillus", f"{FM}/results/lac_vs_uhgg")}
OUT = f"{FM}/results/figS4"


def per_strain_mutualism(resdir, fname):
    """strain -> (n_viable, n_mutualism) at one level."""
    path = f"{resdir}/{fname}.tsv"
    counts = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                ga, gb = float(row["growth_a_alone"]), float(row["growth_b_alone"])
            except (KeyError, ValueError):
                continue
            if not (ga > VIABLE and gb > VIABLE):
                continue
            s = row["species_a"]
            v, m = counts.get(s, (0, 0))
            counts[s] = (v + 1, m + (row.get("interaction_type") == "mutualism"))
    return {s: 100.0 * m / v for s, (v, m) in counts.items() if v}


def boot_mean(x, rng):
    """B bootstrap means of x, resampling elements with replacement."""
    idx = rng.integers(0, len(x), size=(B, len(x)))
    return x[idx].mean(axis=1)


def ci(v):
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    # ── 1. per-strain mutualism at L5 and L6 ──────────────────────────────────
    per = {}   # sys_key -> {strain: {"L5": pct, "L6": pct}}
    for k, (_, resdir) in SYS.items():
        lv = {lab: per_strain_mutualism(resdir, f) for lab, f in LEVELS.items()}
        strains = sorted(set(lv["L5"]) & set(lv["L6"]))
        per[k] = {s: {lab: lv[lab][s] for lab in LEVELS} for s in strains}

    with open(f"{OUT}/per_strain_mutualism_L5_L6.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["system", "strain", "L5_mut_pct", "L6_mut_pct"])
        for k in ("akk", "lac"):
            for s, d in per[k].items():
                w.writerow([k, s, f"{d['L5']:.2f}", f"{d['L6']:.2f}"])

    vals = {k: {lab: np.array([per[k][s][lab] for s in per[k]]) for lab in LEVELS}
            for k in SYS}

    # ── 2. genus mean + bootstrap CI ──────────────────────────────────────────
    rows_ci = []
    for k, (label, _) in SYS.items():
        for lab in ("L5", "L6"):
            x = vals[k][lab]
            lo, hi = ci(boot_mean(x, rng))
            rows_ci.append([label, lab, f"{x.mean():.3f}", f"{lo:.3f}", f"{hi:.3f}", len(x)])
    with open(f"{OUT}/genus_ci_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["system", "level", "mean", "lo", "hi", "n"])
        w.writerows(rows_ci)

    # ── 3. generalist vs specialist at L5 ─────────────────────────────────────
    a5, l5 = vals["akk"]["L5"], vals["lac"]["L5"]
    ba, bl = boot_mean(a5, rng), boot_mean(l5, rng)
    ratio_ci = ci(bl / ba)
    diff_ci = ci(bl - ba)
    mw = mannwhitneyu(l5, a5, alternative="two-sided", method="exact")

    tests = [["Lactobacillus/Akkermansia mutualism ratio at L5",
              f"{l5.mean()/a5.mean():.2f}", f"{ratio_ci[0]:.2f}", f"{ratio_ci[1]:.2f}",
              "percentile bootstrap over strains", "", "",
              len(a5), len(l5), f"B={B}, seed={SEED}"],
             ["Lactobacillus - Akkermansia mutualism difference at L5 (percentage points)",
              f"{l5.mean()-a5.mean():.1f}", f"{diff_ci[0]:.1f}", f"{diff_ci[1]:.1f}",
              "percentile bootstrap over strains", "", "",
              len(a5), len(l5), f"B={B}, seed={SEED}"],
             ["Generalist vs specialist at L5 (rank test across strains)",
              "", "", "", "Mann-Whitney U, two-sided, exact",
              f"{mw.statistic:.0f}", f"{mw.pvalue:.3f}", len(a5), len(l5),
              "marginal: the genera overlap at the strain level"]]

    # ── 4. paired L5 -> L6 decline within each genus ──────────────────────────
    for k, (label, _) in SYS.items():
        d = vals[k]["L6"] - vals[k]["L5"]
        idx = rng.integers(0, len(d), size=(B, len(d)))
        lo, hi = ci(d[idx].mean(axis=1))
        wx = wilcoxon(vals[k]["L5"], vals[k]["L6"], alternative="two-sided", method="exact")
        n_down = int((d < 0).sum())
        tests.append([f"{label}: paired L5 -> L6 change in mutualism (percentage points)",
                      f"{d.mean():.1f}", f"{lo:.1f}", f"{hi:.1f}",
                      "Wilcoxon signed-rank, two-sided, exact",
                      f"{wx.statistic:.0f}", f"{wx.pvalue:.4f}", len(d), len(d),
                      f"{n_down}/{len(d)} strains decline; bootstrap B={B}, seed={SEED}"])

    with open(f"{OUT}/strain_level_tests.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["quantity", "estimate", "ci_lo", "ci_hi", "test", "statistic",
                    "p_value", "n_akkermansia", "n_lactobacillus", "note"])
        w.writerows(tests)

    print(f"wrote {OUT}/per_strain_mutualism_L5_L6.tsv, genus_ci_summary.tsv, strain_level_tests.tsv")
    for t in tests:
        est = f"{t[1]} [{t[2]}, {t[3]}]" if t[1] else ""
        p = f"p = {t[6]}" if t[6] else ""
        print(f"  {t[0][:62]:64} {est:22} {p}")


if __name__ == "__main__":
    main()
