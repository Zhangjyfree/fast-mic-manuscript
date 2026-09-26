#!/usr/bin/env python3
"""Mechanistic control for the Akkermansia glucose crash (Table S10, Discussion).

The Akkermansia models import free glucose and maltose through transporters that
carry no gene association. This script blocks every non-exchange reaction that
moves glucose (cpd00027), maltose (cpd00179) or maltodextrin (cpd11976) from the
extracellular to the cytosolic compartment, checks that Akkermansia monoculture
growth is then identical on L5, L5-glc and L6 (i.e. the added sugars are no longer
usable), and re-screens the blocked models against the UHGG community on the three
media of Figure 4G.

Run from the repository root after extracting test/UHGG/final_gapseq_xml (README).
The fast-mic binary and the media come from the engine checkout, $FASTMIC_ENGINE
(default ../fast-mic).

Output: results/fig4/akk_sugarblock_{L5_pectin,L5glc,L6_resistant_starch}.tsv
        results/fig4/akk_sugarblock_summary.tsv  (viable pairs and mutualism per medium)
Usage : python3 scripts/akk_sugar_block.py [--threads N]
"""
import argparse, csv, glob, os, re, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, default=0)
a = ap.parse_args()

ENGINE = os.environ.get("FASTMIC_ENGINE", "../fast-mic")
BIN = os.path.join(ENGINE, "target", "release", "fast-mic")
MODELS = "test/akk/akk_gapseq_xml"
PARTNERS = "test/UHGG/final_gapseq_xml"
LEVELS = ["L5_pectin", "L5glc", "L6_resistant_starch"]
SUGARS = {"cpd00027": "glucose", "cpd00179": "maltose", "cpd11976": "maltodextrin"}
VIABLE = 1e-4
RXN = re.compile(r'<reaction\b[^>]*\bid="([^"]+)"[^>]*>.*?</reaction>', re.S)

TMP = tempfile.TemporaryDirectory()               # blocked models, removed on exit
blocked_dir = os.path.join(TMP.name, "akk_sugarblock")
os.makedirs(blocked_dir)

# 1. block the free-sugar transporters
for path in sorted(glob.glob(f"{MODELS}/*.xml")):
    x = open(path, encoding="utf-8").read()
    hits = []
    for rid, body in ((m.group(1), m.group(0)) for m in RXN.finditer(x)):
        if rid.startswith(("R_EX_", "R_DM_", "R_SK_")):
            continue
        species = set(re.findall(r'species="M_(cpd\d+)_([ce])0"', body))
        names = [n for c, n in SUGARS.items() if (c, "e") in species and (c, "c") in species]
        if names:
            hits.append((rid, "+".join(names), "gene-supported" if "geneProductAssociation" in body else "no gene support"))
    for rid, _, _ in hits:
        for side in ("lower", "upper"):
            x = re.sub(r'(<reaction\b[^>]*\bid="%s"[^>]*?)fbc:%sFluxBound="[^"]*"' % (re.escape(rid), side),
                       r'\1fbc:%sFluxBound="default_0"' % side, x)
    open(os.path.join(blocked_dir, os.path.basename(path)), "w", encoding="utf-8").write(x)
    print(f"{os.path.basename(path)}: blocked {hits}")

def screen(group2, level, out, threads):
    subprocess.run([BIN, "--group1", blocked_dir, "--group2", group2,
                    "--medium-file", f"{ENGINE}/media/gradient_{level}_gapseq.csv",
                    "--threads", str(threads), "-o", out],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return list(csv.DictReader(open(out), delimiter="\t"))

# 2. the added sugars must no longer change Akkermansia monoculture growth
probe = os.path.join(TMP.name, "probe")
os.makedirs(probe)
first = sorted(f for f in os.listdir(PARTNERS) if f.endswith(".xml"))[0]
os.symlink(os.path.abspath(os.path.join(PARTNERS, first)), os.path.join(probe, first))
growth = {lv: {r["species_a"]: float(r["growth_a_alone"])
               for r in screen(probe, lv, os.path.join(TMP.name, f"probe_{lv}.tsv"), 1)}
          for lv in LEVELS}
for strain in growth["L5_pectin"]:
    g = [growth[lv][strain] for lv in LEVELS]
    if max(g) - min(g) > 1e-6:
        sys.exit(f"{strain}: monoculture growth still changes across L5/L5-glc/L6 {g}")
print("monoculture growth identical on L5, L5-glc and L6 for all strains")

# 3. community screen on the three media
os.makedirs("results/fig4", exist_ok=True)
rows = []
for lv in LEVELS:
    out = f"results/fig4/akk_sugarblock_{lv}.tsv"
    res = screen(PARTNERS, lv, out, a.threads)
    viable = [r for r in res if float(r["growth_a_alone"]) > VIABLE and float(r["growth_b_alone"]) > VIABLE]
    mut = sum(r["interaction_type"] == "mutualism" for r in viable)
    rows.append([lv, len(viable), mut, f"{100 * mut / len(viable):.2f}"])
    print(f"-> {out}  mutualism {100 * mut / len(viable):.2f}% of {len(viable)} viable pairs")
with open("results/fig4/akk_sugarblock_summary.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["medium", "viable_pairs", "mutualistic_pairs", "mutualism_pct"])
    w.writerows(rows)
print("-> results/fig4/akk_sugarblock_summary.tsv")
