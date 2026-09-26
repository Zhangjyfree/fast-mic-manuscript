#!/usr/bin/env python3
"""Objective-function sensitivity (Supplementary Figure S9, Table S12).

Screens L. gasseri ATCC 33323 against the bacterial UHGG species representatives
(the 20 archaeal genomes are excluded, using the UHGG metadata) at L5 and L6,
once with the default lexicographic max-min co-culture objective and once with
the fixed-ratio objective (--fixed-ratio).

Run from the repository root after extracting test/UHGG/final_gapseq_xml (README).
The fast-mic binary and the media come from the engine checkout, $FASTMIC_ENGINE
(default ../fast-mic).

Output: results/figS9/{L5_pectin,L6_resistant_starch}_{lexicographic,fixed_ratio}.tsv
Usage : python3 scripts/figS9_objective_sensitivity.py [--threads N]
"""
import argparse, csv, gzip, os, subprocess, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, default=0)
a = ap.parse_args()
TMP = tempfile.TemporaryDirectory()          # symlink farm, removed when the script exits

ENGINE = os.environ.get("FASTMIC_ENGINE", "../fast-mic")
BIN = os.path.join(ENGINE, "target", "release", "fast-mic")
PROBIOTIC = "test/lac/lac_genomes_faa_gapseq_wdm_xml/L_gasseri_ATCC33323_GCF_000014425.1.xml"
UHGG = "test/UHGG/final_gapseq_xml"
META = "test/UHGG/genomes_metadata_with_gtdb.tsv.gz"

with gzip.open(META, "rt") as fh:
    archaea = {r["Genome"] for r in csv.DictReader(fh, delimiter="\t")
               if r.get("Lineage", "").startswith("d__Archaea")}

# fast-mic takes directories, so link the selected models into the work directory
g1, g2 = os.path.join(TMP.name, "probiotic"), os.path.join(TMP.name, "partners")
for d in (g1, g2):
    os.makedirs(d)
os.symlink(os.path.abspath(PROBIOTIC), os.path.join(g1, os.path.basename(PROBIOTIC)))
n = 0
for f in sorted(os.listdir(UHGG)):
    if f.endswith(".xml") and f[:-4] not in archaea:
        os.symlink(os.path.abspath(os.path.join(UHGG, f)), os.path.join(g2, f)); n += 1
print(f"partners: {n} bacterial UHGG models ({len(archaea & {f[:-4] for f in os.listdir(UHGG)})} archaea excluded)")

os.makedirs("results/figS9", exist_ok=True)
for level in ("L5_pectin", "L6_resistant_starch"):
    for obj, extra in (("lexicographic", []), ("fixed_ratio", ["--fixed-ratio"])):
        out = f"results/figS9/{level}_{obj}.tsv"
        subprocess.run([BIN, "--group1", g1, "--group2", g2,
                        "--medium-file", f"{ENGINE}/media/gradient_{level}_gapseq.csv",
                        "--threads", str(a.threads), "-o", out] + extra,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"-> {out}")
