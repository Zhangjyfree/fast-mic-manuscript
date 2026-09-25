#!/usr/bin/env python
"""A. Scaling timing — fast-mic vs SMETANA on pairwise interaction typing (Fig 1 headline).

Times both tools on the SAME set of n = {100, 500, 1000} probiotic×UHGG pairs, single-
threaded, 3 replicates, on one machine (report the M3 Pro), then extrapolates SMETANA's
per-pair cost to full catalogue scale to quantify "SMETANA infeasible / fast-mic feasible".

Design (apples-to-apples, identical inputs):
  - Unit of work = one species pair -> interaction result.
    SMETANA computes MIP+MRO per pair (global mode, exactly as run_smetana_batch.py).
    fast-mic computes full six-type classification + realized fluxes per pair (does MORE
    per pair, yet is far faster — state this in the legend).
  - Pair set = ALL probiotic strains × a sampled UHGG subset, so the cross-product is
    exactly n pairs (n_uhgg = n / n_prob). Nested: the 50- and 10-UHGG sets are prefixes
    of the 100-UHGG set. Same pairs across replicates (reps measure timing noise only).
  - Single-thread by default (--threads 1) to match the single-threaded COBRApy benchmark
    in Fig 1B. SMETANA is single-threaded per pair regardless.

Run in the conda env `smetana` (needs reframed + smetana + CPLEX/Gurobi):
    python scripts/scaling_timing_smetana.py --sys lac --sizes 100,500,1000 --reps 3
Output:
    results/smetana_comparison/scaling_timing_{sys}_{level}.tsv   (raw per-run rows)
    results/smetana_comparison/scaling_extrapolation_{sys}_{level}.txt  (headline numbers)
"""
import sys, os, csv, time, random, argparse, subprocess, tempfile, statistics, warnings, resource
warnings.filterwarnings("ignore")

FM = "/Users/jingyi/fast-mic"
BIN = f"{FM}/target/release/fast-mic"

PROBDIR = {
    "lac": f"{FM}/test/lac/lac_genomes_faa_gapseq_wdm_xml",
    "akk": f"{FM}/test/akk/akk_gapseq_xml",   # same models as the gradient screen
}
UHGGDIR = f"{FM}/test/UHGG/final_gapseq_xml"
MEDCSV = {
    "L5_pectin": "gradient_L5_pectin_gapseq.csv",
    "L6_resistant_starch": "gradient_L6_resistant_starch_gapseq.csv",
}

ap = argparse.ArgumentParser()
ap.add_argument("--sys", default="lac", choices=["lac", "akk"])
ap.add_argument("--level", default="L6_resistant_starch", choices=list(MEDCSV))
ap.add_argument("--sizes", default="100,500,1000", help="comma-separated pair counts")
ap.add_argument("--reps", type=int, default=3)
ap.add_argument("--threads", type=int, default=1, help="fast-mic threads (1 = fair vs single-threaded SMETANA)")
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--min-growth", type=float, default=0.01)
ap.add_argument("--max-uptake", type=float, default=10.0)
ap.add_argument("--solver", default="cplex", choices=["cplex", "gurobi", "scip"],
                help="reframed MILP backend for SMETANA. reframed auto-picks gurobi>cplex>scip "
                     "if unset, so set explicitly. Note: reframed has NO native HiGHS backend "
                     "(only a slow pulp_highs wrapper, never auto-selected) — SMETANA's performant "
                     "path needs a commercial MILP solver, whereas fast-mic uses free native HiGHS.")
ap.add_argument("--out-suffix", default="",
                help="appended to the output filenames, so a test run cannot overwrite existing results")
ap.add_argument("--smetana-cap-pairs", type=int, default=1000,
                help="skip SMETANA above this pair count (it is too slow); still time fast-mic")
args = ap.parse_args()

SIZES = [int(x) for x in args.sizes.split(",")]
probdir = PROBDIR[args.sys]
med_path = f"{FM}/media/{MEDCSV[args.level]}"
med = [(r["compounds"], float(r["maxFlux"])) for r in csv.DictReader(open(med_path))]

prob_ids = sorted(os.path.splitext(f)[0] for f in os.listdir(probdir) if f.endswith(".xml"))
n_prob = len(prob_ids)
uhgg_ids = sorted(os.path.splitext(f)[0] for f in os.listdir(UHGGDIR) if f.endswith(".xml"))
random.seed(args.seed)
random.shuffle(uhgg_ids)

outdir = f"{FM}/results/smetana_comparison"
os.makedirs(outdir, exist_ok=True)
raw_path = f"{outdir}/scaling_timing_{args.sys}_{args.level}{args.out_suffix}.tsv"
raw = csv.writer(open(raw_path, "w"), delimiter="\t")
raw.writerow(["tool", "n_pairs", "n_prob", "n_uhgg", "threads", "solver", "rep", "wall_s", "sec_per_pair", "peak_rss_mb"])

# ── peak-memory helpers ───────────────────────────────────────────────────────
# ru_maxrss is bytes on macOS/BSD and kilobytes on Linux.
_RSS_SCALE = 1.0 if sys.platform == "darwin" else 1024.0

def rss_mb(ru_maxrss):
    return ru_maxrss * _RSS_SCALE / (1024.0 ** 2)

def self_peak_mb():
    return rss_mb(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

# ── SMETANA (imported lazily so fast-mic-only runs don't require the env) ──────
def smetana_pairs(pairs):
    """Time SMETANA global MIP+MRO over `pairs`; mirror run_smetana_batch.py exactly."""
    from math import inf
    from reframed import load_cbmodel, Environment, set_default_solver
    from reframed.core.model import ReactionType
    from smetana.legacy import Community
    from smetana.smetana import mip_score, mro_score

    set_default_solver(args.solver)   # force CPLEX (else reframed auto-picks gurobi>cplex>scip)
    cache = {}
    def prep(path, mid):
        if mid in cache: return cache[mid]
        m = load_cbmodel(path, flavor="fbc2"); m.id = mid
        m.set_objective({"R_bio1": 1.0}); m.biomass_reaction = "R_bio1"
        if "R_EX_cpd11416_c0" in m.reactions:
            m.reactions["R_EX_cpd11416_c0"].reaction_type = ReactionType.SINK
        cache[mid] = m; return m

    base_peak = self_peak_mb()
    t0 = time.time()
    for a, b in pairs:
        m1 = prep(f"{probdir}/{a}.xml", a); m2 = prep(f"{UHGGDIR}/{b}.xml", b)
        comm = Community("c", [m1, m2], copy_models=True)
        env = Environment()
        for cpd, mf in med:
            env[f"R_EX_M_{cpd}_e0_pool"] = (-mf, inf)
        try:
            mip_score(comm, environment=env, min_growth=args.min_growth, max_uptake=args.max_uptake, verbose=False)
            mro_score(comm, environment=env, min_growth=args.min_growth, max_uptake=args.max_uptake, verbose=False)
        except Exception:
            pass  # failed pairs still cost wall-clock; that is the honest cost
    return time.time() - t0, max(self_peak_mb(), base_peak)

# ── fast-mic (subprocess; times model load + all-pairs typing, single process) ─
def fastmic_pairs(uhgg_subset):
    """Symlink the sampled UHGG models into a temp dir so fast-mic loads ONLY those,
    then run the full prob × subset cross-product (= n pairs). Returns wall seconds."""
    with tempfile.TemporaryDirectory() as td:
        for uid in uhgg_subset:
            os.symlink(f"{UHGGDIR}/{uid}.xml", f"{td}/{uid}.xml")
        out = f"{td}/out.tsv"
        cmd = [BIN, "--group1", probdir, "--group2", td,
               "--medium-file", med_path, "--threads", str(args.threads),
               "-o", out, "--summary"]
        t0 = time.time()
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        _, status, ru = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") else 0
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)
        return time.time() - t0, rss_mb(ru.ru_maxrss)

# ── run ───────────────────────────────────────────────────────────────────────
print(f"system={args.sys} ({n_prob} strains) × UHGG | level={args.level} | threads={args.threads} | reps={args.reps}")
results = {}  # (tool, n) -> [wall_s,...]
mems = {}     # (tool, n) -> [peak_rss_mb,...]
for n in SIZES:
    if n % n_prob != 0:
        print(f"  ! n={n} not divisible by n_prob={n_prob}; rounding n_uhgg")
    n_uhgg = max(1, round(n / n_prob))
    n_actual = n_uhgg * n_prob
    uhgg_subset = uhgg_ids[:n_uhgg]
    pairs = [(p, u) for p in prob_ids for u in uhgg_subset]  # identical set for both tools

    for rep in range(1, args.reps + 1):
        w, mem = fastmic_pairs(uhgg_subset)
        results.setdefault(("fast-mic", n_actual), []).append(w)
        mems.setdefault(("fast-mic", n_actual), []).append(mem)
        raw.writerow(["fast-mic", n_actual, n_prob, n_uhgg, args.threads, "highs", rep,
                      f"{w:.3f}", f"{w/n_actual:.5f}", f"{mem:.1f}"])
        print(f"  [fast-mic] n={n_actual:5d} rep{rep}: {w:8.2f} s  ({w/n_actual*1000:.2f} ms/pair)  peak {mem:.0f} MB")

    if n_actual <= args.smetana_cap_pairs:
        for rep in range(1, args.reps + 1):
            w, mem = smetana_pairs(pairs)
            results.setdefault(("SMETANA", n_actual), []).append(w)
            mems.setdefault(("SMETANA", n_actual), []).append(mem)
            raw.writerow(["SMETANA", n_actual, n_prob, n_uhgg, 1, args.solver, rep,
                          f"{w:.3f}", f"{w/n_actual:.5f}", f"{mem:.1f}"])
            print(f"  [SMETANA ] n={n_actual:5d} rep{rep}: {w:8.2f} s  ({w/n_actual:.3f} s/pair)  peak {mem:.0f} MB")
    else:
        print(f"  [SMETANA ] n={n_actual} > cap {args.smetana_cap_pairs}: skipped (extrapolated instead)")

# ── summary + extrapolation ─────────────────────────────────────────────────────
def med_mm(xs): return (statistics.median(xs), min(xs), max(xs))

lines = []
lines.append(f"# Scaling timing — {args.sys} × UHGG, {args.level}, threads={args.threads}, reps={args.reps}")
lines.append(f"# solvers: fast-mic=HiGHS (native, free) | SMETANA={args.solver.upper()} (reframed has no native HiGHS)")
lines.append(f"# machine: <fill in — Apple M3 Pro>\n")
lines.append(f"{'tool':10s} {'n_pairs':>8s} {'median_s':>10s} {'min_s':>9s} {'max_s':>9s} {'s/pair(median)':>15s}")
sp_per_pair = {}  # tool -> median s/pair at the largest measured n
for (tool, n), xs in sorted(results.items(), key=lambda kv: (kv[0][0], kv[0][1])):
    m, lo, hi = med_mm(xs)
    lines.append(f"{tool:10s} {n:8d} {m:10.2f} {lo:9.2f} {hi:9.2f} {m/n:15.5f}")
    sp_per_pair[tool] = m / n  # last (largest n) wins for each tool

lines.append("")
lines.append(f"{'tool':10s} {'n_pairs':>8s} {'peak_MB(median)':>16s} {'min':>8s} {'max':>8s}")
for (tool, n), xs in sorted(mems.items(), key=lambda kv: (kv[0][0], kv[0][1])):
    m, lo, hi = med_mm(xs)
    lines.append(f"{tool:10s} {n:8d} {m:16.1f} {lo:8.1f} {hi:8.1f}")
lines.append("# fast-mic memory = peak RSS of the fast-mic child process (single process, all pairs).")
lines.append("# SMETANA memory = peak RSS of the Python process running reframed/SMETANA over the same pairs.")

# full-scale extrapolation
n_uhgg_full = len(uhgg_ids)
n_media = 10
study_pairs = n_prob * n_uhgg_full * n_media          # this study's screen (one probiotic panel × UHGG × L0–L9)
catalogue_pairs = n_uhgg_full * (n_uhgg_full - 1) // 2  # full UHGG all-vs-all (single medium)

def humantime(sec):
    for unit, s in [("years", 3.156e7), ("days", 86400.0), ("hours", 3600.0), ("min", 60.0)]:
        if sec >= s: return f"{sec/s:.1f} {unit}"
    return f"{sec:.1f} s"

lines.append("")
lines.append("## Extrapolation (median s/pair × full pair count)")
for label, npair in [(f"this study ({args.sys} {n_prob}×UHGG {n_uhgg_full}×{n_media} media)", study_pairs),
                     (f"full UHGG all-vs-all ({n_uhgg_full} genomes, 1 medium)", catalogue_pairs)]:
    lines.append(f"  {label}: {npair:,} pairs")
    for tool in ("SMETANA", "fast-mic"):
        if tool in sp_per_pair:
            lines.append(f"      {tool:9s}: {humantime(sp_per_pair[tool]*npair)}  "
                         f"({sp_per_pair[tool]:.5f} s/pair)")
    lines.append("")
if "SMETANA" in sp_per_pair and "fast-mic" in sp_per_pair:
    lines.append(f"## Single-thread speedup (s/pair): {sp_per_pair['SMETANA']/sp_per_pair['fast-mic']:.0f}×")

txt = "\n".join(lines)
open(f"{outdir}/scaling_extrapolation_{args.sys}_{args.level}{args.out_suffix}.txt", "w").write(txt + "\n")
print("\n" + txt)
print(f"\nraw -> {raw_path}")
