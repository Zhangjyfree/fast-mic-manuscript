#!/usr/bin/env python
"""Extract per-strain data for Fig 2 v2 (phylogeny + metabolic repertoire + niche breadth).

For each probiotic system, writes two tidy TSVs into results/fig2/:
  fig2_traits_{sys}.tsv  — one row per strain, GEM-reconstruction descriptors:
      reactions, metabolites, genes, transporters, exchanges, pathways, noGPR_reactions
  fig2_growth_{sys}.tsv  — long format strain × prebiotic level × monoculture growth rate
      (growth_a_alone, de-duplicated per strain per level from the gradient co-culture TSVs;
       a strain's monoculture rate is identical across all its pairs at a level)

Trait definitions (namespace-agnostic SBML parse; no external deps):
  reactions   = <reaction> count (all)
  metabolites = <species> count
  genes       = <geneProduct> definitions (fbc)
  transporters= non-exchange reactions whose species span >=2 compartments (c0/e0/p0)
  exchanges   = reactions with id R_EX_/R_DM_/R_SK_
  pathways    = <group> count (gapseq subsystems)
  noGPR       = non-exchange reactions lacking a geneProductAssociation (gapfill+spontaneous proxy)

Usage:
  python scripts/extract_fig2_traits.py            # both lac and akk
  python scripts/extract_fig2_traits.py lac        # one system
"""
import sys, os, csv, glob
import xml.etree.ElementTree as ET

FM = "."   # repo root
# Models used for the gradient screen: Lactobacillus gap-filled on Western diet +
# mucin, Akkermansia gap-filled on the Akkermansia minimal medium (Methods).
SYS_DIRS = {
    "lac": (f"{FM}/test/lac/lac_genomes_faa_gapseq_wdm_xml", f"{FM}/results/lac_vs_uhgg"),
    "akk": (f"{FM}/test/akk/akk_gapseq_xml", f"{FM}/results/akk_vs_uhgg"),
}
LEVELS = [("L0_base","L0"),("L1_inulin","L1"),("L2_fos","L2"),("L3_gos","L3"),
          ("L4_xos","L4"),("L5_pectin","L5"),("L6_resistant_starch","L6"),
          ("L7_bglucan","L7"),("L8_hmo","L8"),("L9_mos","L9")]

def lname(tag): return tag.split('}')[-1]

def gem_stats(path):
    root = ET.parse(path).getroot()
    genes = sum(1 for e in root.iter() if lname(e.tag) == 'geneProduct')
    pathways = sum(1 for e in root.iter() if lname(e.tag) == 'group')
    mets = sum(1 for e in root.iter() if lname(e.tag) == 'species')
    rxn = exch = transport = noGPR = 0
    for r in root.iter():
        if lname(r.tag) != 'reaction':
            continue
        rxn += 1
        rid = r.get('id', '')
        comps, has_gpa = set(), False
        for sub in r.iter():
            lt = lname(sub.tag)
            if lt == 'geneProductAssociation':
                has_gpa = True
            elif lt == 'speciesReference':
                sp = sub.get('species', '')
                for c in ('_e0', '_c0', '_p0'):
                    if sp.endswith(c):
                        comps.add(c)
        if rid.startswith(('R_EX_', 'R_DM_', 'R_SK_')):
            exch += 1
            continue
        if len(comps) >= 2:
            transport += 1
        if not has_gpa:
            noGPR += 1
    return dict(reactions=rxn, metabolites=mets, genes=genes, transporters=transport,
                exchanges=exch, pathways=pathways, noGPR_reactions=noGPR)

def growth_by_level(resdir):
    """strain -> {level_label: growth_a_alone}."""
    out = {}
    for fname, lab in LEVELS:
        p = f"{resdir}/{fname}.tsv"
        if not os.path.exists(p):
            continue
        seen = {}
        for row in csv.DictReader(open(p), delimiter="\t"):
            s = row["species_a"]
            if s not in seen:
                try:
                    seen[s] = float(row["growth_a_alone"])
                except (ValueError, KeyError):
                    seen[s] = float("nan")
        for s, g in seen.items():
            out.setdefault(s, {})[lab] = g
    return out

def run(sys_key):
    gdir, resdir = SYS_DIRS[sys_key]
    os.makedirs(f"{FM}/results/fig2", exist_ok=True)
    strains = sorted(os.path.splitext(os.path.basename(f))[0] for f in glob.glob(f"{gdir}/*.xml"))

    # traits
    tpath = f"{FM}/results/fig2/fig2_traits_{sys_key}.tsv"
    cols = ["reactions","metabolites","genes","transporters","exchanges","pathways","noGPR_reactions"]
    with open(tpath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["strain"] + cols)
        for s in strains:
            st = gem_stats(f"{gdir}/{s}.xml")
            w.writerow([s] + [st[c] for c in cols])
    print(f"[{sys_key}] {len(strains)} strains -> {tpath}")

    # growth (long)
    gmap = growth_by_level(resdir)
    gpath = f"{FM}/results/fig2/fig2_growth_{sys_key}.tsv"
    levels_present = [lab for (_, lab) in LEVELS]
    with open(gpath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["strain", "level", "growth"])
        for s in strains:
            for lab in levels_present:
                g = gmap.get(s, {}).get(lab, "")
                if g != "":
                    w.writerow([s, lab, f"{g:.6f}"])
    n_lev = len(set(l for d in gmap.values() for l in d))
    print(f"[{sys_key}] growth over {n_lev} levels -> {gpath}")

if __name__ == "__main__":
    keys = sys.argv[1:] or ["lac", "akk"]
    for k in keys:
        run(k)
