#!/usr/bin/env python3
"""
Benchmark single-species FBA / pFBA using COBRApy.
Replicates fast-mic's algorithm exactly:
  1. Apply medium with tiered uptake limits + cofactor preservation (medium.rs)
  2. FBA: maximise biomass
  3. CycleFreeFlux + pFBA (single LP): fix exchanges, cap biomass, minimise Σ|v_i|

Usage:
    python benchmark_cobra.py <media_db.tsv> <medium_name> <model1.xml> [...]
    python benchmark_cobra.py <media_db.tsv> <medium_name> --model-list <file.txt>
    python benchmark_cobra.py <media_db.tsv> <medium_name> --model-list f.txt --workers 4
"""

import sys
import csv
import glob
import time
import resource
import gc
import os
import signal
from pathlib import Path
from multiprocessing import Pool, TimeoutError as MPTimeoutError

try:
    import cobra
except ImportError:
    print("ERROR: COBRApy not installed.  pip install cobra", file=sys.stderr)
    sys.exit(1)

# ── Use HiGHS solver to match fast-mic (Rust uses HiGHS via good_lp) ──
#
# COBRApy 0.31 / optlang 1.9 reach HiGHS through `optlang.hybrid_interface`,
# which requires BOTH `highspy` (HiGHS Python bindings) AND `osqp` (the
# hybrid interface keeps OSQP as a QP fallback even when only LP is used).
#
# Falls back to GLPK with a clear hint if either package is missing,
# so the benchmark still runs — just on a slower LP backend.
def _configure_solver():
    try:
        from optlang import hybrid_interface  # noqa: F401 — needs highspy + osqp
        cobra.Configuration().solver = hybrid_interface
        return "HiGHS (via optlang.hybrid_interface)"
    except Exception as e:
        print(
            f"WARNING: HiGHS solver not available ({e}); falling back to GLPK.\n"
            "         To match fast-mic's LP backend, install both:\n"
            "             pip install highspy osqp",
            file=sys.stderr,
        )
        try:
            cobra.Configuration().solver = "glpk"
        except Exception:
            pass
        return f"GLPK (fallback) — {cobra.Configuration().solver}"

_SOLVER_NAME = _configure_solver()

# ── Constants matching src/cobra.rs ──
MIN_VIABLE_GROWTH = 1e-4
NUMERICAL_TOL     = 1e-5

# Per-medium wall-clock timeout (seconds). GLPK can hang on degenerate LPs.
# Override with the MODEL_TIMEOUT_S environment variable.
MODEL_TIMEOUT_S = int(os.environ.get("MODEL_TIMEOUT_S", "120"))

# ── Tiered uptake limits matching src/medium.rs MediumBounds::default() ──
CARBON_SOURCE_LIMIT = 10.0
AMINO_ACID_LIMIT    = 1.0
NUCLEOBASE_LIMIT    = 0.5
COFACTOR_LIMIT      = 0.1

# ── Compound classification lists (ported from src/medium.rs) ──
AMINO_ACIDS_BIGG = {
    "ala__l","arg__l","asn__l","asp__l","cys__l","gln__l","glu__l","gly",
    "his__l","ile__l","leu__l","lys__l","met__l","phe__l","pro__l","ser__l",
    "thr__l","trp__l","tyr__l","val__l","ala__d","arg__d","asn__d","asp__d",
    "cys__d","gln__d","glu__d","his__d","ile__d","leu__d","lys__d","met__d",
    "phe__d","pro__d","ser__d","thr__d","trp__d","tyr__d","val__d",
    "orn","citr","hcys__l","gaba","taurine","cgly","5mthf",
}
AMINO_ACIDS_SEED = {
    "cpd00023","cpd00033","cpd00035","cpd00039","cpd00041","cpd00051",
    "cpd00053","cpd00054","cpd00060","cpd00065","cpd00066","cpd00069",
    "cpd00084","cpd00107","cpd00117","cpd00119","cpd00129","cpd00132",
    "cpd00135","cpd00156","cpd00161","cpd00186","cpd00268","cpd00320",
    "cpd00322","cpd00345","cpd00549","cpd00550","cpd00567","cpd00586",
    "cpd00587","cpd00637","cpd01017","cpd19016",
}
NUCLEOBASES_BIGG = {
    "ura","gua","ade","cyt","thy","hyxn","uri","guo","ado","cyd","thd","ino",
    "xan","xanthine","amp","gmp","cmp","ump","imp","adn","gsn","csn","urd",
    "duri","dguo","dado","dcyt","dtmp","dump",
}
NUCLEOBASES_SEED = {
    "cpd00018","cpd00046","cpd00091","cpd00092","cpd00114","cpd00126",
    "cpd00128","cpd00182","cpd00207","cpd00249","cpd00298","cpd00299",
    "cpd00307","cpd00309","cpd00311","cpd00412","cpd00654",
}
COFACTORS_BIGG = {
    "fol","thf","5mthf","mlthf","nmn","nad","nadp","fmn","fad","ribflv",
    "cbl1","cbl2","adocbl","mecbl","btn","pnto__r","pnto__s","4ppcys","pap",
    "pyam5p","pydx","pydxn","pydx5p","thm","thmpp","pheme","sheme","hem",
    "mqn7","mqn8","mql8","q8","q8h2","lipoate","lipoamide","spmd","ptrc","cadav",
}
COFACTORS_SEED = {
    "cpd00003","cpd00006","cpd00015","cpd00016","cpd00028","cpd00045",
    "cpd00050","cpd00056","cpd00087","cpd00104","cpd00118","cpd00125",
    "cpd00166","cpd00215","cpd00220","cpd00263","cpd00264","cpd00305",
    "cpd00355","cpd00393","cpd00423","cpd00493","cpd00541","cpd00557",
    "cpd00635","cpd00644","cpd02666","cpd11606","cpd15499","cpd15500",
    "cpd15561",
}


def _met_base(met_id: str) -> str:
    s = met_id.lower()
    s = s[2:] if s.startswith("m_") else s
    for suffix in ("_e0","_e","_c0","_c","_p0","_p"):
        if s.endswith(suffix):
            s = s[:-len(suffix)]
            break
    return s


def uptake_limit_for_met(met_id: str) -> float:
    base = _met_base(met_id)
    if base in AMINO_ACIDS_BIGG or base in AMINO_ACIDS_SEED:
        return AMINO_ACID_LIMIT
    if base in NUCLEOBASES_BIGG or base in NUCLEOBASES_SEED:
        return NUCLEOBASE_LIMIT
    if base in COFACTORS_BIGG or base in COFACTORS_SEED:
        return COFACTOR_LIMIT
    return CARBON_SOURCE_LIMIT


def is_cofactor_met(met_id: str) -> bool:
    base = _met_base(met_id)
    return base in COFACTORS_BIGG or base in COFACTORS_SEED


# ============================================================
# Media DB
# ============================================================

def parse_media_db(path):
    db = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            name, desc, compound = row[0].strip(), row[1].strip(), row[2].strip()
            if name not in db:
                db[name] = (desc, set())
            db[name][1].add(compound)
    return db


def parse_medium_csv(path):
    """Parse a gapseq SEED-format medium CSV (compounds,name,maxFlux).
    Returns (compounds_set, per_compound_max_flux_dict).
    """
    compounds, per_bounds = set(), {}
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader, None)  # header
        for row in reader:
            if len(row) < 3:
                continue
            cpd = row[0].strip()
            if not cpd:
                continue
            try:
                mx = float(row[2].strip())
            except ValueError:
                continue
            if mx > 0:
                compounds.add(cpd)
                per_bounds[cpd] = mx
    return compounds, per_bounds


def expand_medium_compounds(base_compounds):
    medium = set()
    for raw in base_compounds:
        medium.add(raw)
        wp = raw[2:] if raw.startswith("M_") else raw
        medium.add(wp)
        medium.add(f"M_{wp}")
        root = wp
        for suffix in ("_e","_e0","_c","_c0","_p","_p0"):
            if root.endswith(suffix):
                root = root[:-len(suffix)]
                break
        medium.add(root)
        medium.add(f"M_{root}")
        for ext in ("_e","_e0"):
            medium.add(f"{root}{ext}")
            medium.add(f"M_{root}{ext}")
    return medium


def met_in_medium(met_id, expanded):
    if met_id in expanded:
        return True
    s = met_id[2:] if met_id.startswith("M_") else met_id
    if s in expanded:
        return True
    root = s
    for suffix in ("_e","_e0"):
        if root.endswith(suffix):
            root = root[:-len(suffix)]
            break
    else:
        return False
    return root in expanded or f"M_{root}" in expanded


# ============================================================
# Exchange detection
# ============================================================

EXCHANGE_EXCLUDES = ("demand","DM_","biosynthesis","transcription",
                     "replication","SN_","SK_","sink")

def find_exchanges(model):
    try:
        native = list(model.exchanges)
        if len(native) >= 5:
            return native
    except Exception:
        pass
    prefix = [r for r in model.reactions
              if r.id.startswith("EX_") or r.id.startswith("R_EX_")]
    if prefix:
        return prefix
    return [r for r in model.reactions
            if len(r.metabolites) == 1
            and not any(ex in r.id for ex in EXCHANGE_EXCLUDES)]


# ============================================================
# Medium application
# ============================================================

def apply_medium_fastmic(model, base_compounds, per_compound_bounds=None):
    expanded = expand_medium_compounds(base_compounds)
    exchanges = find_exchanges(model)
    n_opened = n_closed = 0
    for rxn in exchanges:
        mets = list(rxn.metabolites.keys())
        if not mets:
            continue
        met     = mets[0]
        met_id  = met.id
        stoich  = rxn.metabolites[met]
        in_med  = met_in_medium(met_id, expanded)
        # Per-compound CSV override → exact maxFlux; otherwise tiered classification
        if per_compound_bounds:
            base_id = _met_base(met_id)
            limit = per_compound_bounds.get(base_id, uptake_limit_for_met(met_id))
        else:
            limit = uptake_limit_for_met(met_id)
        is_cof  = is_cofactor_met(met_id)
        sbml_lb = rxn.lower_bound
        sbml_ub = rxn.upper_bound
        if stoich < 0:
            if in_med:
                rxn.lower_bound = -limit;  n_opened += 1
            elif is_cof and sbml_lb < -1e-9:
                n_opened += 1
            else:
                rxn.lower_bound = 0.0;  n_closed += 1
        else:
            if in_med:
                rxn.upper_bound = limit;  n_opened += 1
            elif is_cof and sbml_ub > 1e-9:
                n_opened += 1
            else:
                rxn.upper_bound = 0.0;  n_closed += 1
    return exchanges, n_opened, n_closed


# ============================================================
# CycleFreeFlux + pFBA
# ============================================================

def cycle_free_pfba(model, biomass_rxn, exchanges):
    exchange_ids = {r.id for r in exchanges}
    with model:
        model.objective = biomass_rxn
        model.objective.direction = "max"
        sol = model.optimize()
    if sol.status != "optimal":
        return 0.0
    fba_growth = sol.fluxes.get(biomass_rxn.id, 0.0)
    if fba_growth < MIN_VIABLE_GROWTH:
        return 0.0
    with model:
        for rxn in exchanges:
            v = sol.fluxes.get(rxn.id, 0.0)
            rxn.lower_bound = v - NUMERICAL_TOL
            rxn.upper_bound = v + NUMERICAL_TOL
        biomass_rxn.lower_bound = 0.0
        biomass_rxn.upper_bound = fba_growth
        internal = [r for r in model.reactions
                    if r.id not in exchange_ids and r.id != biomass_rxn.id]
        iface = model.solver.interface
        abs_vars = {}
        for rxn in internal:
            av = iface.Variable(f"__abs_{rxn.id}", lb=0.0)
            model.solver.add(av)
            abs_vars[rxn.id] = av
            rv = model.solver.variables[rxn.id]
            model.solver.add(iface.Constraint(av - rv, lb=0.0,
                                               name=f"__c1_{rxn.id}"))
            model.solver.add(iface.Constraint(av + rv, lb=0.0,
                                               name=f"__c2_{rxn.id}"))
        model.solver.objective = iface.Objective(
            sum(abs_vars.values()), direction="min", name="cff_obj")
        status = model.optimize().status
        if status != "optimal":
            return fba_growth
        return model.solver.variables[biomass_rxn.id].primal


# ============================================================
# Worker (runs in a subprocess — safe to kill on timeout)
# ============================================================

def _worker(args):
    """Run one model under every medium. Returns a LIST of result dicts,
    one per (model, medium). The model is loaded once; each medium keeps its
    own growth rate (NOT averaged), so the correctness check compares fast-mic
    vs COBRApy per medium. A single medium yields a length-1 list with the bare
    model id (back-compatible)."""
    path, media = args
    model_id_base = Path(path).stem
    multi = len(media) > 1

    if hasattr(signal, "SIGALRM"):
        def _alarm(sig, frame):
            raise TimeoutError(f"SIGALRM after {MODEL_TIMEOUT_S}s")
        signal.signal(signal.SIGALRM, _alarm)

    try:
        rows = []
        for (base_compounds, per_bounds, label) in media:
            # Re-read the model for each medium: COBRApy's CycleFreeFlux mutates
            # the solver (adds __abs_* vars), so a fresh model per medium is the
            # robust choice. fast-mic loads once on its side; the growth values
            # (what the correctness check compares) are identical either way.
            if hasattr(signal, "SIGALRM"):
                signal.alarm(MODEL_TIMEOUT_S)        # per-medium timeout
            t0 = time.perf_counter()
            model = cobra.io.read_sbml_model(path)
            load_s = time.perf_counter() - t0

            n_mets  = len(model.metabolites)
            n_rxns  = len(model.reactions)
            n_genes = len(model.genes)
            biomass_rxn = None
            for rxn in model.reactions:
                rid = rxn.id.lower()
                if "biomass" in rid or "biomass" in (rxn.name or "").lower() \
                        or rid == "growth":
                    biomass_rxn = rxn
                    break
            biomass_id = biomass_rxn.id if biomass_rxn else "NONE"

            exchanges, n_opened, n_closed = apply_medium_fastmic(model, base_compounds, per_bounds)
            t1 = time.perf_counter()
            growth = cycle_free_pfba(model, biomass_rxn, exchanges) if biomass_rxn else 0.0
            pfba_s = time.perf_counter() - t1
            mid = f"{model_id_base}__{label}" if multi else model_id_base
            rows.append({
                "model_id": mid, "n_mets": n_mets, "n_rxns": n_rxns,
                "n_genes": n_genes, "biomass_rxn": biomass_id,
                "growth": growth, "load_s": load_s, "pfba_s": pfba_s,
                "n_opened": n_opened, "n_closed": n_closed,
                "timed_out": False, "error": None,
            })

        if hasattr(signal, "SIGALRM"):
            signal.alarm(0)
        return rows

    except Exception as exc:
        if hasattr(signal, "SIGALRM"):
            signal.alarm(0)
        return [{
            "model_id": model_id_base, "n_mets": 0, "n_rxns": 0,
            "n_genes": 0, "biomass_rxn": "NONE",
            "growth": 0.0, "load_s": 0.0, "pfba_s": 0.0,
            "n_opened": 0, "n_closed": 0,
            "timed_out": isinstance(exc, TimeoutError),
            "error": str(exc),
        }]


# ============================================================
# Memory helpers
# ============================================================

def get_peak_memory_mb():
    u = resource.getrusage(resource.RUSAGE_SELF)
    return u.ru_maxrss / (1024*1024) if sys.platform == "darwin" \
           else u.ru_maxrss / 1024


# ============================================================
# Main
# ============================================================

def main():
    # ── Parse args ──
    #
    # Two invocation modes:
    #   (A) <media_db.tsv> <medium_name>  [models or --model-list]
    #   (B) --medium-file <CSV>           [models or --model-list]
    #
    # Mode B uses a gapseq-style SEED CSV (compounds,name,maxFlux) and
    # is the right choice for ModelSEED / gapseq SBMLs.
    raw_argv = sys.argv[1:]

    workers      = 1
    medium_file  = None
    media_list   = None
    raw_paths    = []
    leftover     = []

    # Two-pass parse: extract flagged args first so positionals are unambiguous.
    i = 0
    while i < len(raw_argv):
        a = raw_argv[i]
        if a == "--model-list":
            i += 1
            with open(raw_argv[i]) as fh:
                for line in fh:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        raw_paths.append(line)
        elif a in ("--workers", "--threads"):
            i += 1
            workers = int(raw_argv[i])
            if workers <= 0:
                workers = os.cpu_count() or 1
        elif a == "--medium-file":
            i += 1
            medium_file = raw_argv[i]
        elif a == "--media-list":
            i += 1
            media_list = raw_argv[i]
        else:
            leftover.append(a)
        i += 1

    if medium_file is not None or media_list is not None:
        media_db_path = None
        medium_name   = None
        pos_models    = leftover
    else:
        if len(leftover) < 2:
            print(
                "Usage: benchmark_cobra.py <media_db.tsv> <medium_name> "
                "<model1.xml> [...]\n"
                "       benchmark_cobra.py <media_db.tsv> <medium_name> "
                "--model-list <file.txt> [--workers N]\n"
                "       benchmark_cobra.py --medium-file <medium.csv> "
                "--model-list <file.txt> [--workers N]",
                file=sys.stderr,
            )
            sys.exit(1)
        media_db_path = leftover[0]
        medium_name   = leftover[1]
        pos_models    = leftover[2:]

    # Expand any glob patterns in positional model args
    for p in pos_models:
        if "*" in p or "?" in p:
            raw_paths.extend(sorted(glob.glob(p)))
        else:
            raw_paths.append(p)

    model_paths = raw_paths
    if not model_paths:
        print("ERROR: no model paths", file=sys.stderr)
        sys.exit(1)

    print(f"Read {len(model_paths)} model paths", file=sys.stderr)

    # ── Resolve medium / media ──
    # `media` is a list of (base_compounds, per_bounds, label). Each model is
    # loaded once and evaluated under every medium (per-model inner loop).
    # For a single medium this is a length-1 list (original behaviour).
    if media_list is not None:
        with open(media_list) as fh:
            mpaths = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        if not mpaths:
            print(f"ERROR: media-list '{media_list}' is empty", file=sys.stderr); sys.exit(1)
        media = []
        for mp in mpaths:
            bc, pb = parse_medium_csv(mp)
            media.append((bc, pb, Path(mp).stem))
        per_bounds = media[0][1]
        medium_label = f"{len(media)} media (per-model loop): " + ", ".join(m[2] for m in media)
    elif medium_file is not None:
        base_compounds, per_bounds = parse_medium_csv(medium_file)
        media = [(base_compounds, per_bounds, Path(medium_file).stem)]
        medium_label = Path(medium_file).stem
    else:
        media_db = parse_media_db(media_db_path)
        if medium_name not in media_db:
            print(f"Medium '{medium_name}' not found. Available: {sorted(media_db)}",
                  file=sys.stderr)
            sys.exit(1)
        desc, base_compounds = media_db[medium_name]
        media = [(base_compounds, None, medium_name)]
        per_bounds = None
        medium_label = f"{medium_name} ({desc})"

    print(f"Medium: {medium_label}", file=sys.stderr)
    print(f"Models: {len(model_paths)}", file=sys.stderr)
    print(f"COBRApy: {cobra.__version__}  solver: {cobra.Configuration().solver}",
          file=sys.stderr)
    print(f"Algorithm: FBA → CycleFreeFlux+pFBA (fast-mic style)", file=sys.stderr)
    print(f"Uptake limits: C={CARBON_SOURCE_LIMIT} aa={AMINO_ACID_LIMIT} "
          f"nuc={NUCLEOBASE_LIMIT} cof={COFACTOR_LIMIT}"
          f"{'  (per-compound overrides from CSV)' if per_bounds else ''}",
          file=sys.stderr)
    print(f"Workers: {workers}  Per-model timeout: {MODEL_TIMEOUT_S}s",
          file=sys.stderr)
    print(file=sys.stderr)

    total_start = time.perf_counter()
    work = [(p, media) for p in model_paths]
    results_ordered = [None] * len(model_paths)
    n_timeout = 0
    n_error = 0

    if workers == 1:
        # Single-process path — simpler, easier to debug
        for idx, item in enumerate(work):
            path = item[0]
            model_id = Path(path).stem
            print(f"  [{idx+1}/{len(model_paths)}] {model_id}",
                  file=sys.stderr, flush=True)
            rows = _worker(item)
            r0 = rows[0]
            if r0["timed_out"]:
                n_timeout += 1
                print(f"    TIMEOUT (>{MODEL_TIMEOUT_S}s) — skipped",
                      file=sys.stderr)
            elif r0["error"]:
                n_error += 1
                print(f"    ERROR: {r0['error']}", file=sys.stderr)
            else:
                gmin = min(x['growth'] for x in rows); gmax = max(x['growth'] for x in rows)
                print(f"    {len(rows)} media  mets={r0['n_mets']} rxns={r0['n_rxns']} "
                      f"biomass={r0['biomass_rxn']} growth={gmin:.4f}..{gmax:.4f}  "
                      f"load={r0['load_s']:.3f}s  pfba={sum(x['pfba_s'] for x in rows):.3f}s",
                      file=sys.stderr)
            results_ordered[idx] = rows
    else:
        # Multi-process path with per-model timeout
        with Pool(processes=workers) as pool:
            futures = [(idx, pool.apply_async(_worker, (item,)))
                       for idx, item in enumerate(work)]
            for idx, fut in futures:
                path = work[idx][0]
                model_id = Path(path).stem
                print(f"  [{idx+1}/{len(model_paths)}] {model_id}",
                      file=sys.stderr, flush=True)
                try:
                    rows = fut.get(timeout=MODEL_TIMEOUT_S + 5)
                except MPTimeoutError:
                    rows = [{
                        "model_id": model_id, "n_mets": 0, "n_rxns": 0,
                        "n_genes": 0, "biomass_rxn": "NONE",
                        "growth": 0.0, "load_s": 0.0, "pfba_s": 0.0,
                        "n_opened": 0, "n_closed": 0,
                        "timed_out": True, "error": "pool timeout",
                    }]
                r0 = rows[0]
                if r0["timed_out"]:
                    n_timeout += 1
                    print(f"    TIMEOUT (>{MODEL_TIMEOUT_S}s) — skipped",
                          file=sys.stderr)
                elif r0["error"]:
                    n_error += 1
                    print(f"    ERROR: {r0['error']}", file=sys.stderr)
                else:
                    gmin = min(x['growth'] for x in rows); gmax = max(x['growth'] for x in rows)
                    print(f"    {len(rows)} media  mets={r0['n_mets']} rxns={r0['n_rxns']} "
                          f"biomass={r0['biomass_rxn']} growth={gmin:.4f}..{gmax:.4f}  "
                          f"load={r0['load_s']:.3f}s  pfba={sum(x['pfba_s'] for x in rows):.3f}s",
                          file=sys.stderr)
                results_ordered[idx] = rows

    total_s     = time.perf_counter() - total_start
    peak_mem_mb = get_peak_memory_mb()
    # Flatten: each model contributed a list of per-medium rows.
    valid = [r for rows in results_ordered if rows
             for r in rows if not r["timed_out"] and not r["error"]]
    n_valid = len(valid)

    # ── TSV to stdout (one row per (model, medium)) ──
    print("model_id\tn_metabolites\tn_reactions\tn_genes\t"
          "biomass_rxn\tgrowth_rate\tload_time_s\tpfba_time_s")
    for rows in results_ordered:
        if not rows:
            continue
        for r in rows:
            print(f"{r['model_id']}\t{r['n_mets']}\t{r['n_rxns']}\t{r['n_genes']}\t"
                  f"{r['biomass_rxn']}\t{r['growth']:.6f}\t"
                  f"{r['load_s']:.6f}\t{r['pfba_s']:.6f}")

    # ── Summary to stderr ──
    tl = sum(r["load_s"] for r in valid)
    tp = sum(r["pfba_s"] for r in valid)
    n  = n_valid or 1
    print(file=sys.stderr)
    print("=" * 44, file=sys.stderr)
    print("  COBRApy CycleFreeFlux+pFBA Benchmark", file=sys.stderr)
    print("=" * 44, file=sys.stderr)
    print(f"Models:          {len(model_paths)}", file=sys.stderr)
    print(f"Completed:       {n_valid}", file=sys.stderr)
    if n_timeout:
        print(f"Timed out:       {n_timeout}  (>{MODEL_TIMEOUT_S}s each)",
              file=sys.stderr)
    if n_error:
        print(f"Errors:          {n_error}", file=sys.stderr)
    if n_valid:
        print(f"Total load:      {tl:.1f} s  ({tl/n:.3f} s/model)",
              file=sys.stderr)
        print(f"Total CFF+pFBA:  {tp:.1f} s  ({tp/n:.3f} s/model)",
              file=sys.stderr)
    print(f"Total wall time: {total_s:.1f} s", file=sys.stderr)
    print(f"Peak memory:     {peak_mem_mb:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
