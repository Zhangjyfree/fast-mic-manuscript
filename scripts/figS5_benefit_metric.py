#!/usr/bin/env python3
"""Relative versus absolute benefit metric (Supplementary Figure S5).

For each system and gradient level, mutualism is scored two ways over the viable
pairs (both monoculture growth rates > 1e-4 h^-1):
  relative  both benefits beta = (mu_co - mu_alone) / mu_alone > 1e-3 (as in the paper)
  absolute  both mu_co - mu_alone > 1e-3 h^-1 (no mu_alone denominator)
together with the mean monoculture growth and the mean absolute co-culture benefit.

Input : results/{akk,lac}_vs_uhgg/L*.tsv
Output: results/figS5/absolute_vs_relative.tsv
Usage : python3 scripts/figS5_benefit_metric.py        (from the repository root)
"""
import csv, os, statistics as st

LEVELS = [("L0", "L0_base"), ("L1", "L1_inulin"), ("L2", "L2_fos"), ("L3", "L3_gos"),
          ("L4", "L4_xos"), ("L5", "L5_pectin"), ("L6", "L6_resistant_starch"),
          ("L7", "L7_bglucan"), ("L8", "L8_hmo"), ("L9", "L9_mos")]
EPS_B = 1e-3; EPS_ABS = 1e-3; VIA = 1e-4
OUT = "results/figS5/absolute_vs_relative.tsv"

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["system", "level", "viable", "mu_alone_mean", "mutualism_relative_pct",
                "mutualism_absolute_pct", "mean_abs_benefit"])
    for system, d in [("Akkermansia", "results/akk_vs_uhgg"), ("Lactobacillus", "results/lac_vs_uhgg")]:
        for lv, fn in LEVELS:
            via = 0; mu = []; rel = 0; ab = 0; de = []
            for r in csv.DictReader(open(os.path.join(d, fn + ".tsv")), delimiter="\t"):
                ga, gb = float(r["growth_a_alone"]), float(r["growth_b_alone"])
                if ga <= VIA or gb <= VIA:
                    continue
                via += 1
                da, db = float(r["growth_a_co"]) - ga, float(r["growth_b_co"]) - gb
                mu.append((ga + gb) / 2); de.append((da + db) / 2)
                if float(r["benefit_a"]) > EPS_B and float(r["benefit_b"]) > EPS_B:
                    rel += 1
                if da > EPS_ABS and db > EPS_ABS:
                    ab += 1
            w.writerow([system, lv, via, f"{st.mean(mu):.4f}", f"{100*rel/via:.2f}",
                        f"{100*ab/via:.2f}", f"{st.mean(de):.4f}"])
print(f"-> {OUT}")
