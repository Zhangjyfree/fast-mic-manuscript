#!/usr/bin/env python
"""Fig 5 enrichment: condition-resolved cross-feeding prevalence + metabolite class.

For each system, across L0-L9, among MUTUALISTIC pairs, computes how often each
probiotic-exported metabolite (a_to_b, probiotic = species_a) appears.
Outputs results/fig5/fig5_prevalence_{sys}.tsv with columns:
  cpd, name, class, level, n_mut_pairs, n_export, prevalence(%)

Metabolite name from ModelSEED compounds.tsv; class by keyword on the name
(organic acid / nucleoside / sugar / amino acid / alcohol-diol / other).
"""
import sys, os, csv, gzip, glob
csv.field_size_limit(1 << 24)

FM = "."   # repo root
SYS_RES = {"lac": f"{FM}/results/lac_vs_uhgg", "akk": f"{FM}/results/akk_vs_uhgg"}
LEVELS = [("L0_base","L0"),("L1_inulin","L1"),("L2_fos","L2"),("L3_gos","L3"),
          ("L4_xos","L4"),("L5_pectin","L5"),("L6_resistant_starch","L6"),
          ("L7_bglucan","L7"),("L8_hmo","L8"),("L9_mos","L9")]

def load_cpd_names():
    m = {}
    with open(f"{FM}/compounds.tsv") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            m[row["id"]] = row["name"]
    return m

AMINO = {"alanine","arginine","asparagine","aspartate","cysteine","glutamate","glutamine",
         "glycine","histidine","isoleucine","leucine","lysine","methionine","phenylalanine",
         "proline","serine","threonine","tryptophan","tyrosine","valine"}
def classify(name):
    n = name.lower()
    if any(k in n for k in ("inosine","xanthosine","guanosine","cytidine","adenosine","uridine",
                            "thymidine","nucleoside","xanthine","hypoxanthine")):
        return "nucleoside"
    if any(k in n for k in ("lactate","malate","succinate","fumarate","formate","propionate",
                            "acetate","pyruvate","citrate","oxaloacetate","2-oxoglutarate",
                            "butyrate","valerate","glutarate","glycolate","glyoxylate")):
        return "organic acid"
    if any(k in n for k in ("ribose","xylose","arabinose","glucose","fructose","mannose",
                            "galactose","rhamnose","fucose","pentose","xylulose","ribulose",
                            "maltose","cellobiose")):
        return "sugar"
    if any(k in n for k in ("ethanol","methanol","propanediol","glycerol","butanediol","butanol",
                            "propanol","mannitol","sorbitol")):
        return "alcohol/diol"
    if any(a in n for a in AMINO):
        return "amino acid"
    return "other"

def run(sys_key, names):
    resdir = SYS_RES[sys_key]
    os.makedirs(f"{FM}/results/fig5", exist_ok=True)
    out = f"{FM}/results/fig5/fig5_prevalence_{sys_key}.tsv"
    w = csv.writer(open(out, "w"), delimiter="\t")
    w.writerow(["cpd","name","class","level","n_mut_pairs","n_export","prevalence"])
    for fname, lab in LEVELS:
        gzp, plainp = f"{resdir}/{fname}.full.tsv.gz", f"{resdir}/{fname}.full.tsv"
        if os.path.exists(gzp):
            opener = lambda: gzip.open(gzp, "rt")
        elif os.path.exists(plainp):
            opener = lambda: open(plainp)
        else:
            continue
        n_mut = 0
        counts = {}
        with opener() as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                if row.get("interaction_type") != "mutualism":
                    continue
                n_mut += 1
                mets = row.get("a_to_b_metabolites", "") or ""
                seen = set()
                for tok in mets.split(";"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    cpd = tok.replace("M_", "").rsplit("_e0", 1)[0].rsplit("_c0", 1)[0]
                    seen.add(cpd)
                for cpd in seen:
                    counts[cpd] = counts.get(cpd, 0) + 1
        for cpd, c in counts.items():
            nm = names.get(cpd, cpd)
            w.writerow([cpd, nm, classify(nm), lab, n_mut, c,
                        f"{100.0*c/n_mut:.2f}" if n_mut else "0"])
    print(f"[{sys_key}] -> {out}")

if __name__ == "__main__":
    names = load_cpd_names()
    for k in (sys.argv[1:] or ["lac","akk"]):
        run(k, names)
