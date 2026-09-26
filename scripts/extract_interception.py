#!/usr/bin/env python
"""Higher-order 'interception' analysis — how much could a third party erode a pairwise cross-feed?

Answers Reviewer 1.2 (pairwise only / higher-order effects) WITHOUT needing a 3-member
growth objective, which is not well posed: for N>=3 the Pareto set of growth allocations
is high-dimensional and, absent abundance data, no allocation rule is principled; the
six-type classification is also intrinsically pairwise.

Instead we bound the effect empirically: for every metabolite that fast-mic predicts to be
cross-fed in mutualistic pairs, we count how many OTHER community members carry an uptake
(exchange) reaction for it. A metabolite that most of the community can import is liable to
be intercepted before it reaches the intended partner (pairwise mutualism = an upper bound);
a metabolite few members can import is robust to higher-order competition.

Output: results/figS7/interception.tsv
  cpd, name, class, system, prevalence (mean % of mutualistic pairs, from fig5),
  n_uptakers, interception_frac (= n_uptakers / n_community_models)

Usage: python3 scripts/extract_interception.py
"""
import re, glob, os, csv
from collections import defaultdict

FM = "."
UHGG = f"{FM}/test/UHGG/final_gapseq_xml"
OUT = f"{FM}/results/figS7/interception.tsv"
EX = re.compile(rb'id="R_EX_(cpd\d+)_e0"')

def community_uptake_counts():
    """cpd -> number of community models carrying an exchange reaction for it."""
    files = sorted(glob.glob(f"{UHGG}/*.xml"))
    counts = defaultdict(int)
    for f in files:
        with open(f, "rb") as fh:
            for cpd in set(EX.findall(fh.read())):
                counts[cpd.decode()] += 1
    return counts, len(files)

def main():
    counts, n_models = community_uptake_counts()
    print(f"scanned {n_models} community models; {len(counts)} distinct exchangeable compounds")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    for sys_key, label in (("akk", "Akkermansia"), ("lac", "Lactobacillus")):
        p = f"{FM}/results/fig5/fig5_prevalence_{sys_key}.tsv"
        if not os.path.exists(p):
            print(f"  missing {p} — run extract_fig5_crossfeed.py first"); continue
        agg = defaultdict(list); meta = {}
        for r in csv.DictReader(open(p), delimiter="\t"):
            agg[r["cpd"]].append(float(r["prevalence"]))
            meta[r["cpd"]] = (r["name"], r["class"])
        for cpd, prevs in agg.items():
            name, cls = meta[cpd]
            n_up = counts.get(cpd, 0)
            rows.append([cpd, name, cls, label, f"{sum(prevs)/len(prevs):.2f}",
                         n_up, f"{n_up/n_models:.4f}"])
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cpd", "name", "class", "system", "prevalence", "n_uptakers", "interception_frac"])
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {OUT}")

if __name__ == "__main__":
    main()
