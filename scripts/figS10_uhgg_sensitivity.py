#!/usr/bin/env python3
"""Sensitivity of the gradient pattern to the UHGG genome set (Supplementary Figure S10).

Mutualism (% of viable pairs) per system and level over three partner sets: all
species representatives, high-quality genomes only (CheckM completeness >= 90% and
contamination <= 5%, UHGG metadata), and all genomes except Patescibacteriota (CPR).

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
CPR = ("Patescibacteriota", "Patescibacteria")

meta = {}
with gzip.open(META, "rt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        try:
            comp, cont = float(row["Completeness"]), float(row["Contamination"])
        except ValueError:
            comp, cont = 0, 100
        phy = next((p[3:] for p in row.get("Lineage", "").split(";") if p.startswith("p__")), "")
        meta[row["Genome"]] = (comp, cont, phy)

def mut_pct(rows, keep):
    v = m = 0
    for r in rows:
        if float(r["growth_a_alone"]) <= VIA or float(r["growth_b_alone"]) <= VIA:
            continue
        if not keep(r["species_b"]):
            continue
        v += 1; m += r["interaction_type"] == "mutualism"
    return (100 * m / v, v) if v else (float("nan"), 0)

def hq(b):
    c = meta.get(b); return bool(c) and c[0] >= 90 and c[1] <= 5

def not_cpr(b):
    c = meta.get(b); return bool(c) and c[2] not in CPR

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["system", "level", "full_pct", "full_n", "hq_pct", "hq_n", "noCPR_pct", "noCPR_n"])
    for sysn, d in [("Akkermansia", "results/akk_vs_uhgg"), ("Lactobacillus", "results/lac_vs_uhgg")]:
        for lv, fn in LEVELS:
            rows = list(csv.DictReader(open(f"{d}/{fn}.tsv"), delimiter="\t"))
            fp, fn_ = mut_pct(rows, lambda b: True); hp, hn = mut_pct(rows, hq); cp, cn = mut_pct(rows, not_cpr)
            w.writerow([sysn, lv, f"{fp:.2f}", fn_, f"{hp:.2f}", hn, f"{cp:.2f}", cn])
print(f"-> {OUT}")
