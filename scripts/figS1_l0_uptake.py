#!/usr/bin/env python3
"""Substrate basis of Lactobacillus viability at L0 (Supplementary Figure S1, Table S3).

For each Lactobacillus-group model, the monoculture CycleFreeFlux + pFBA solution
on the L0 base medium is taken from fast-mic (--target-reaction on every exchange
reaction) and uptake fluxes are summed by substrate class:
  Amino acids          the 20 proteinogenic L-amino acids
  Mucin amino sugars   N-acetylglucosamine, glucosamine, N-acetylneuraminate,
                       N-acetylmannosamine, L-fucose
  Other carbon         every other carbon-containing substrate (CO2 excluded)

Run from the repository root after extracting test/UHGG/final_gapseq_xml (README;
one partner model is needed only because fast-mic runs pairwise). The fast-mic
binary, the L0 medium and ModelSEED compounds.tsv come from the engine checkout,
$FASTMIC_ENGINE (default ../fast-mic).

Output: results/figS1/lac_L0_uptake.tsv   (strain, growth [h^-1], category, uptake_flux)
Usage : python3 scripts/figS1_l0_uptake.py
"""
import csv, glob, os, re, subprocess, tempfile

ENGINE = os.environ.get("FASTMIC_ENGINE", "../fast-mic")
BIN = os.path.join(ENGINE, "target", "release", "fast-mic")
COMPOUNDS = os.path.join(ENGINE, "media", "compounds.tsv")
MODELS = "test/lac/lac_genomes_faa_gapseq_wdm_xml"
PARTNER = "test/UHGG/final_gapseq_xml/MGYG000000001.xml"
MEDIUM = os.path.join(ENGINE, "media", "gradient_L0_base_gapseq.csv")
OUT = "results/figS1/lac_L0_uptake.tsv"

AMINO = {"cpd00035", "cpd00051", "cpd00132", "cpd00041", "cpd00084", "cpd00023", "cpd00053",
         "cpd00033", "cpd00119", "cpd00322", "cpd00107", "cpd00039", "cpd00060", "cpd00066",
         "cpd00129", "cpd00054", "cpd00161", "cpd00065", "cpd00069", "cpd00156"}
MUCIN = {"cpd00122", "cpd00276", "cpd00232", "cpd00492", "cpd00751"}
INORGANIC_C = {"cpd00011", "cpd00242"}          # CO2, bicarbonate
CATS = ["Amino acids", "Other carbon", "Mucin amino sugars"]

formula = {r["id"]: r["formula"] for r in csv.DictReader(open(COMPOUNDS), delimiter="\t")}
def has_carbon(cpd):
    return re.search(r"C(?![a-z])", formula.get(cpd, "") or "") is not None

def category(cpd):
    if cpd in AMINO: return "Amino acids"
    if cpd in MUCIN: return "Mucin amino sugars"
    if cpd not in INORGANIC_C and has_carbon(cpd): return "Other carbon"
    return None

TMP = tempfile.TemporaryDirectory()          # removed when the script exits
WORK = TMP.name
g1, g2 = f"{WORK}/probiotic", f"{WORK}/partner"
os.makedirs(g1); os.makedirs(g2)
os.symlink(os.path.abspath(PARTNER), f"{g2}/{os.path.basename(PARTNER)}")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
rows = []
for path in sorted(glob.glob(f"{MODELS}/*.xml")):
    strain = re.sub(r"_GCF_.*$", "", os.path.basename(path)[:-4])
    for f in os.listdir(g1): os.remove(f"{g1}/{f}")
    os.symlink(os.path.abspath(path), f"{g1}/{os.path.basename(path)}")
    ex = re.findall(r'<reaction\b[^>]*\bid="(R_EX_cpd\d+_e0)"', open(path).read())
    cmd = [BIN, "--group1", g1, "--group2", g2, "--medium-file", MEDIUM, "--threads", "1",
           "-o", f"{WORK}/out.tsv"] + [a for r in ex for a in ("--target-reaction", r)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    r = next(csv.DictReader(open(f"{WORK}/out.tsv"), delimiter="\t"))
    tot = dict.fromkeys(CATS, 0.0)
    for rid in ex:
        v = float(r[f"{rid}__alone_a"])
        cat = category(rid[5:-3])
        if v < -1e-9 and cat:
            tot[cat] += -v
    for cat in CATS:
        rows.append([strain, f"{float(r['growth_a_alone']):.6f}", cat, round(tot[cat], 3)])
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["strain", "growth", "category", "uptake_flux"]); w.writerows(rows)
print(f"-> {OUT}")
