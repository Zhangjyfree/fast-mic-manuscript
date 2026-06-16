#!/usr/bin/env python3
"""Regenerate results/figures_paper/crossfeed_landscape_table.tsv from the per-pair
*.full.tsv outputs.

For each system (gut only: Akkermansia x UHGG, Lactobacillus x UHGG) this pools all
mutualistic pairs across the L0-L9 gradient and, for every cross-fed metabolite,
reports the number of mutualistic pairs that exchange it in EITHER direction
(a->b union b->a, de-duplicated per pair). This matches the prevalence figures
used in the manuscript (e.g. Akkermansia methanol 49%, succinate 30%).

Inputs : results/<system>/<level>.full.tsv  (columns a_to_b_metabolites,
         b_to_a_metabolites, interaction_type, growth_a_alone, growth_b_alone)
         compounds.tsv  (cpd id -> name)
Output : results/figures_paper/crossfeed_landscape_table.tsv
         columns: sys_label, metabolite, cpd, n_pairs, n_mut, frac

Usage  : python3 scripts/make_crossfeed_table.py
"""
import csv
import os
from collections import defaultdict

THRESH = 1e-4
LEVELS = ["L0_base", "L1_inulin", "L2_fos", "L3_gos", "L4_xos", "L5_pectin",
          "L6_resistant_starch", "L7_bglucan", "L8_hmo", "L9_mos"]
SYSTEMS = [("akk_vs_uhgg", "Akkermansia × UHGG (gut)"),
           ("lac_vs_uhgg", "Lactobacillus × UHGG (gut)")]
OUT = "results/figures_paper/crossfeed_landscape_table.tsv"


def strip(tok):
    return tok.replace("M_", "").replace("_e0", "")


def cpd_names():
    name = {}
    with open("compounds.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("id") and row.get("name"):
                name[row["id"]] = row["name"]
    return name


def tally(sysdir):
    """Return (n_mut, {cpd: n_pairs}) pooled over all levels (both directions)."""
    n_mut = 0
    cnt = defaultdict(int)
    for lv in LEVELS:
        path = f"results/{sysdir}/{lv}.full.tsv"
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if (float(r["growth_a_alone"]) > THRESH and
                        float(r["growth_b_alone"]) > THRESH and
                        r["interaction_type"] == "mutualism"):
                    n_mut += 1
                    mets = set()
                    for col in ("a_to_b_metabolites", "b_to_a_metabolites"):
                        s = r.get(col, "")
                        if s and s.strip():
                            mets.update(strip(x) for x in s.split(";") if x)
                    for m in mets:
                        cnt[m] += 1
    return n_mut, cnt


def main():
    names = cpd_names()
    rows = []
    for sysdir, label in SYSTEMS:
        n_mut, cnt = tally(sysdir)
        for cpd, n in sorted(cnt.items(), key=lambda kv: kv[1], reverse=True):
            rows.append({
                "sys_label": label,
                "metabolite": names.get(cpd, cpd),
                "cpd": cpd,
                "n_pairs": n,
                "n_mut": n_mut,
                "frac": n / n_mut if n_mut else 0.0,
            })
        print(f"{label}: n_mut={n_mut}, {len(cnt)} metabolites")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sys_label", "metabolite", "cpd",
                                          "n_pairs", "n_mut", "frac"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
