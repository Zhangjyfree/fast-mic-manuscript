#!/usr/bin/env python3
"""Gap-fill-candidate knockout (Supplementary Figure S12).

Tests whether predicted cooperation depends on gap-filled metabolism. For each
probiotic model, a reaction is a gap-fill candidate when it

  * has no gene association (no fbc:geneProductAssociation),
  * is not an exchange / demand / sink reaction and not the biomass reaction,
  * is not a transport reaction (species in >= 2 compartments; ModelSEED leaves
    transporters GPR-free irrespective of gap-filling, so removing them would
    abolish cross-feeding mechanically rather than test gap-filling), and
  * carries zero flux (|v| < 1e-9) in the monoculture parsimonious FBA solution
    on the test medium, computed with COBRApy (dormant in monoculture, so it
    could only be recruited for cross-feeding).

Candidates are blocked (both flux bounds set to 0), monoculture growth is checked
to be unchanged, the knockout models are re-screened against the community on
the same medium, and each control mutualistic pair is followed to its new class.

Requires COBRApy and the fast-mic engine checkout ($FASTMIC_ENGINE or --engine,
default ../fast-mic), which provides the binary and the media.

Usage (from the repository root):
  python3 scripts/gapfill_knockout.py --system Akkermansia \
      --models test/akk/akk_gapseq_xml \
      --partners test/UHGG/final_gapseq_xml \
      --control results/akk_vs_uhgg/L5_pectin.tsv \
      --medium L5_pectin --threads 12 --append

  and the same with --system Lactobacillus --models test/lac/lac_genomes_faa_gapseq_wdm_xml
  --control results/lac_vs_uhgg/L5_pectin.tsv. --append replaces that system's row in
  results/figS12/R15_summary.tsv.
"""
import argparse, csv, os, subprocess, sys, tempfile

VIABLE = 1e-4          # viability threshold (h^-1), as in the main analysis
ZERO_FLUX = 1e-9       # |v| below this counts as zero in the monoculture pFBA solution
SUMMARY_COLS = ["system", "ctrl_mut_pct", "aprime_mut_pct", "n_ctrl_mut", "retained_mut",
                "to_competition", "to_neutral_commensal", "retention_pct"]


def run_engine(engine, g1, g2, medium, out, threads):
    cmd = [f"{engine}/target/release/fast-mic", "--group1", g1, "--group2", g2,
           "--medium-file", f"{engine}/media/gradient_{medium}_gapseq.csv",
           "--threads", str(threads), "-o", out]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=os.environ.get("FASTMIC_ENGINE", "../fast-mic"),
                    help="fast-mic engine checkout (default: $FASTMIC_ENGINE or ../fast-mic)")
    ap.add_argument("--system", required=True, help="label written to the summary, e.g. Akkermansia")
    ap.add_argument("--models", required=True, help="directory of probiotic SBML models")
    ap.add_argument("--partners", required=True, help="directory of community SBML models")
    ap.add_argument("--control", required=True, help="control per-pair TSV on the same medium")
    ap.add_argument("--medium", default="L5_pectin")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--workdir", default=None,
                    help="scratch directory (default: a temporary directory, removed afterwards)")
    ap.add_argument("--summary", default="results/figS12/R15_summary.tsv")
    ap.add_argument("--append", action="store_true", help="replace/add this system's row in --summary")
    a = ap.parse_args()

    engine = a.engine
    tmp = None if a.workdir else tempfile.TemporaryDirectory()
    wd = a.workdir or tmp.name
    ko_dir = f"{wd}/ko_models"
    os.makedirs(ko_dir, exist_ok=True)

    # 1. identify and block gap-fill candidates model by model (COBRApy pFBA)
    import warnings
    warnings.filterwarnings("ignore")
    import cobra
    from cobra.flux_analysis import pfba
    med = {r["compounds"]: float(r["maxFlux"])
           for r in csv.DictReader(open(f"{engine}/media/gradient_{a.medium}_gapseq.csv"))}
    for fn in sorted(f for f in os.listdir(a.models) if f.endswith(".xml")):
        m = cobra.io.read_sbml_model(f"{a.models}/{fn}")
        exch = {r.id for r in m.reactions if r.id.startswith("EX_")}
        m.medium = {f"EX_{c}_e0": v for c, v in med.items() if f"EX_{c}_e0" in exch}
        mu0 = m.slim_optimize()
        flux = pfba(m).fluxes
        n = 0
        for r in m.reactions:
            rid = r.id
            if rid.startswith(("EX_", "DM_", "SK_")) or "bio" in rid.lower():
                continue                                   # boundary and biomass reactions
            if r.gene_reaction_rule.strip():
                continue                                   # has gene evidence
            if abs(flux.get(rid, 0.0)) >= ZERO_FLUX:
                continue                                   # active in monoculture
            if len({met.compartment for met in r.metabolites}) > 1:
                continue                                   # transporter: keep
            r.bounds = (0.0, 0.0); n += 1
        mu1 = m.slim_optimize()
        cobra.io.write_sbml_model(m, f"{ko_dir}/{fn}")
        print(f"{fn}: {n} gap-fill candidates blocked; monoculture growth {mu0:.6f} -> {mu1:.6f}")
        if abs(mu1 - mu0) > 1e-6:
            sys.exit(f"monoculture growth changed for {fn}; aborting")

    # 2. re-screen the knockout models against the community
    ko_rows = run_engine(engine, ko_dir, a.partners, a.medium, f"{wd}/ko_{a.medium}.tsv", a.threads)

    # 3. compare with the control screen
    def viable(r):
        return float(r["growth_a_alone"]) > VIABLE and float(r["growth_b_alone"]) > VIABLE
    ctrl = {(r["species_a"], r["species_b"]): r["interaction_type"]
            for r in csv.DictReader(open(a.control), delimiter="\t") if viable(r)}
    ko = {(r["species_a"], r["species_b"]): r["interaction_type"] for r in ko_rows if viable(r)}
    ctrl_mut = [k for k, t in ctrl.items() if t == "mutualism"]
    fate = [ko.get(k, "nonviable") for k in ctrl_mut]
    retained = sum(t == "mutualism" for t in fate)
    to_comp = sum(t == "competition" for t in fate)
    out = dict(system=a.system,
               ctrl_mut_pct=f"{100 * len(ctrl_mut) / len(ctrl):.2f}",
               aprime_mut_pct=f"{100 * sum(t == 'mutualism' for t in ko.values()) / len(ko):.2f}",
               n_ctrl_mut=len(ctrl_mut), retained_mut=retained, to_competition=to_comp,
               to_neutral_commensal=len(ctrl_mut) - retained - to_comp,
               retention_pct=f"{100 * retained / len(ctrl_mut):.1f}")
    print("\t".join(SUMMARY_COLS)); print("\t".join(str(out[c]) for c in SUMMARY_COLS))

    if a.append:
        rows = []
        if os.path.exists(a.summary):
            rows = [r for r in csv.DictReader(open(a.summary), delimiter="\t") if r["system"] != a.system]
        rows.append(out)
        rows.sort(key=lambda r: r["system"])
        with open(a.summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print(f"-> {a.summary}")


if __name__ == "__main__":
    main()
