#!/usr/bin/env python3
"""Gene support of reactions versus cross-feeding flux (Supplementary Figure S11).

Per probiotic panel:
  Intracellular reactions  % of non-boundary, single-compartment reactions with a
                           gene association (fbc:geneProductAssociation)
  Transport reactions      the same for reactions whose species span >= 2
                           compartments (the transporter definition of Figure 2)
  Cross-feeding flux       mean gene-supported fraction of cross-feeding flux over
                           all viable pairs, pooled across L0-L9
Model percentages are averaged over the strains of a panel (each strain weighted
equally). Exchange, demand, sink and biomass reactions are excluded.

Input : test/akk/akk_gapseq_xml, test/lac/lac_genomes_faa_gapseq_wdm_xml,
        results/{akk,lac}_vs_uhgg/L*.tsv
Output: results/figS11/gpr_coverage.tsv
Usage : python3 scripts/figS11_gpr_coverage.py        (from the repository root)
"""
import csv, glob, os, re, statistics as st

PANELS = [("Akkermansia", "test/akk/akk_gapseq_xml", "results/akk_vs_uhgg"),
          ("Lactobacillus", "test/lac/lac_genomes_faa_gapseq_wdm_xml", "results/lac_vs_uhgg")]
LEVELS = ["L0_base", "L1_inulin", "L2_fos", "L3_gos", "L4_xos", "L5_pectin",
          "L6_resistant_starch", "L7_bglucan", "L8_hmo", "L9_mos"]
VIA = 1e-4
OUT = "results/figS11/gpr_coverage.tsv"
RXN = re.compile(r'<reaction\b[^>]*\bid="([^"]+)"[^>]*>(.*?)</reaction>', re.S)
COMP = re.compile(r'species="M_[^"]*?_([a-z])\d*"')

def model_coverage(path):
    x = open(path, encoding="utf-8").read()
    n = {"intra": [0, 0], "trans": [0, 0]}
    for rid, body in RXN.findall(x):
        if rid.startswith(("R_EX_", "R_DM_", "R_SK_")) or "bio" in rid.lower():
            continue
        k = "trans" if len(set(COMP.findall(body))) > 1 else "intra"
        n[k][0] += "geneProductAssociation" in body
        n[k][1] += 1
    return 100 * n["intra"][0] / n["intra"][1], 100 * n["trans"][0] / n["trans"][1]

def crossfeed_support(resdir):
    vals = []
    for lv in LEVELS:
        for r in csv.DictReader(open(f"{resdir}/{lv}.tsv"), delimiter="\t"):
            if float(r["growth_a_alone"]) > VIA and float(r["growth_b_alone"]) > VIA:
                try:
                    v = float(r["gene_supported_fraction"])
                except ValueError:
                    continue
                if v == v:                      # skip NaN
                    vals.append(v)
    return 100 * st.mean(vals)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["panel", "category", "gene_support_pct"])
    for panel, mdir, rdir in PANELS:
        cov = [model_coverage(f) for f in sorted(glob.glob(f"{mdir}/*.xml"))]
        w.writerow([panel, "Intracellular reactions", f"{st.mean(c[0] for c in cov):.1f}"])
        w.writerow([panel, "Transport reactions", f"{st.mean(c[1] for c in cov):.1f}"])
        w.writerow([panel, "Cross-feeding flux", f"{crossfeed_support(rdir):.1f}"])
print(f"-> {OUT}")
