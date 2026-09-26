#!/usr/bin/env python3
"""Aggregate the loop-removal validation by prebiotic level (Supplementary Table S5).

Input : results/tableS5/cff_biomass_deviation.tsv (bench-cff-deviation output:
        1,000 UHGG models x 10 media; see README)
Output: results/tableS5/tableS5_by_level.tsv
        level, growing evaluations, max and mean |deviation| (1/h) over growing ones
Usage : python3 scripts/tableS5_summary.py            (from the repository root)
"""
import csv

IN, OUT = "results/tableS5/cff_biomass_deviation.tsv", "results/tableS5/tableS5_by_level.tsv"
LEVELS = ["L0_base", "L1_inulin", "L2_fos", "L3_gos", "L4_xos", "L5_pectin",
          "L6_resistant_starch", "L7_bglucan", "L8_hmo", "L9_mos"]

dev = {lv: [] for lv in LEVELS}
for r in csv.DictReader(open(IN), delimiter="\t"):
    if r["viable"] == "true":
        dev[r["medium"].replace("gradient_", "").replace("_gapseq", "")].append(abs(float(r["deviation"])))
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["level", "growing_evaluations", "max_abs_deviation", "mean_abs_deviation"])
    for lv in LEVELS:
        d = dev[lv]
        w.writerow([lv.split("_")[0], len(d), f"{max(d):.3e}", f"{sum(d)/len(d):.3e}"])
    alld = [v for lv in LEVELS for v in dev[lv]]
    w.writerow(["All", len(alld), f"{max(alld):.3e}", f"{sum(alld)/len(alld):.3e}"])
print(f"-> {OUT}")
