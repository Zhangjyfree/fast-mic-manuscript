#!/usr/bin/env python3
"""Sensitivity of the gradient pattern to the UHGG genome set (Supplementary Figure S10).

The 3,238 community models already pass CheckM2 (completeness >= 90%, contamination
<= 5%) and exclude Patescibacteriota (CPR) (Methods). As a stricter quality check,
mutualism (% of viable pairs) per system and level is recomputed over the subset
of partners that also pass the same thresholds under the original CheckM estimates
distributed with the UHGG catalogue. Partners absent from the UHGG metadata
cannot be checked and are left out of the stricter set.

Input : results/{akk,lac}_vs_uhgg/L*.tsv, test/UHGG/genomes_metadata_with_gtdb.tsv.gz
Output: results/figS10/uhgg_sensitivity.tsv
Usage : python3 scripts/figS10_uhgg_sensitivity.py     (from the repository root)
"""
import csv, gzip, os

META = "test/UHGG/genomes_metadata_with_gtdb.tsv.gz"
OUT = "results/figS10/uhgg_sensitivity.tsv"
LEVELS = [("L0", "L0_base"), ("L1", "L1_inulin"), ("L2", "L2_fos"), ("L3", "L3_gos"),
          ("L4", "L4_xos"), ("L5", "L5_pectin"), ("L6", "L6_resistant_starch"),
          ("L7", "L7_bglucan"), ("L8", "L8_hmo"), ("L9", "L9_mos")]
VIA = 1e-4

meta = {}
with gzip.open(META, "rt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        try:
            meta[row["Genome"]] = (float(row["Completeness"]), float(row["Contamination"]))
        except ValueError:
            pass

def mut_pct(rows, keep):
    v = m = 0
    for r in rows:
        if float(r["growth_a_alone"]) <= VIA or float(r["growth_b_alone"]) <= VIA:
            continue
        if not keep(r["species_b"]):
            continue
        v += 1; m += r["interaction_type"] == "mutualism"
    return (100 * m / v, v) if v else (float("nan"), 0)

def strict(b):
    c = meta.get(b); return bool(c) and c[0] >= 90 and c[1] <= 5

os.makedirs(os.path.dirname(OUT), exist_ok=True)
partners = set()
with open(OUT, "w") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["system", "level", "full_pct", "full_n", "strict_pct", "strict_n"])
    for sysn, d in [("Akkermansia", "results/akk_vs_uhgg"), ("Lactobacillus", "results/lac_vs_uhgg")]:
        for lv, fn in LEVELS:
            rows = list(csv.DictReader(open(f"{d}/{fn}.tsv"), delimiter="\t"))
            partners |= {r["species_b"] for r in rows}
            fp, fn_ = mut_pct(rows, lambda b: True); sp, sn = mut_pct(rows, strict)
            w.writerow([sysn, lv, f"{fp:.2f}", fn_, f"{sp:.2f}", sn])
kept = sum(strict(b) for b in partners)
print(f"partners {len(partners)}; stricter set {kept}; "
      f"without UHGG metadata {sum(b not in meta for b in partners)}")
print(f"-> {OUT}")
