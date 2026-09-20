#!/usr/bin/env python
"""Fig 6 enrichment: per-strain bootstrap CI, and per-strain x partner-phylum cooperation.

Outputs (results/fig6/):
  fig6_bootstrap_{sys}.tsv  strain, mean_mut, ci_lo, ci_hi   (mean over L0-L9 of per-level
                            mutualism% among viable pairs; 95% CI by resampling partners)
  fig6_phylum_{sys}.tsv     strain, phylum, mean_mut, mean_n  (mean over levels of mutualism%
                            among viable partners of each GTDB phylum)
Point estimate of mean_mut also feeds the mechanistic scatter (merged with fig2 traits).
"""
import sys, os, csv, gzip
import numpy as np
csv.field_size_limit(1 << 24)

FM = "."   # repo root
SYS_RES = {"lac": f"{FM}/results/lac_vs_uhgg", "akk": f"{FM}/results/akk_vs_uhgg"}
LEVELS = [("L0_base","L0"),("L1_inulin","L1"),("L2_fos","L2"),("L3_gos","L3"),
          ("L4_xos","L4"),("L5_pectin","L5"),("L6_resistant_starch","L6"),
          ("L7_bglucan","L7"),("L8_hmo","L8"),("L9_mos","L9")]
VIABLE = 1e-4
B = 1000
rng = np.random.default_rng(42)

def load_phylum():
    m = {}
    with gzip.open(f"{FM}/test/UHGG/genomes_metadata_with_gtdb.tsv.gz", "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            tax = row.get("GTDB_Classification", "")
            ph = "Unknown"
            for seg in tax.split(";"):
                if seg.startswith("p__"):
                    ph = seg[3:] or "Unknown"; break
            m[row["Genome"]] = ph
    return m

def read_level(sys_key, fname):
    """yield (strain, partner, viable_bool, mut_bool)."""
    for ext in (".tsv",):
        p = f"{SYS_RES[sys_key]}/{fname}{ext}"
        if not os.path.exists(p):
            return
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    ga = float(row["growth_a_alone"]); gb = float(row["growth_b_alone"])
                except (KeyError, ValueError):
                    continue
                viable = ga > VIABLE and gb > VIABLE
                yield (row["species_a"], row["species_b"], viable,
                       row.get("interaction_type") == "mutualism")

def run(sys_key, phylum):
    os.makedirs(f"{FM}/results/fig6", exist_ok=True)
    # gather: strain -> {partner -> [mut per level or nan if not viable]}
    n_lev = len(LEVELS)
    strains = {}
    for li, (fname, lab) in enumerate(LEVELS):
        for strain, partner, viable, mut in read_level(sys_key, fname):
            d = strains.setdefault(strain, {})
            arr = d.get(partner)
            if arr is None:
                arr = [np.nan] * n_lev; d[partner] = arr
            arr[li] = (1.0 if mut else 0.0) if viable else np.nan

    # ── bootstrap CI of mean-over-levels mutualism% ──
    bpath = f"{FM}/results/fig6/fig6_bootstrap_{sys_key}.tsv"
    bw = csv.writer(open(bpath, "w"), delimiter="\t"); bw.writerow(["strain","mean_mut","ci_lo","ci_hi"])
    for strain, pdict in strains.items():
        M = np.array(list(pdict.values()))  # (n_partner, n_lev), 1/0/nan
        with np.errstate(invalid="ignore"):
            per_lvl = np.nanmean(M, axis=0) * 100.0
        point = np.nanmean(per_lvl)
        n_p = M.shape[0]
        boot = np.empty(B)
        for b in range(B):
            idx = rng.integers(0, n_p, n_p)
            with np.errstate(invalid="ignore"):
                pl = np.nanmean(M[idx], axis=0) * 100.0
            boot[b] = np.nanmean(pl)
        lo, hi = np.nanpercentile(boot, [2.5, 97.5])
        bw.writerow([strain, f"{point:.3f}", f"{lo:.3f}", f"{hi:.3f}"])
    print(f"[{sys_key}] bootstrap -> {bpath}")

    # ── per-strain x phylum mean mutualism ──
    ppath = f"{FM}/results/fig6/fig6_phylum_{sys_key}.tsv"
    pw = csv.writer(open(ppath, "w"), delimiter="\t"); pw.writerow(["strain","phylum","mean_mut","mean_n"])
    for strain, pdict in strains.items():
        # per level: group viable partners by phylum
        # accum[phylum] = [sum_mut_per_level, n_viable_per_level] arrays
        acc = {}
        for partner, arr in pdict.items():
            ph = phylum.get(partner, "Unknown")
            a = acc.setdefault(ph, [np.zeros(n_lev), np.zeros(n_lev)])
            v = np.array(arr)
            viable_mask = ~np.isnan(v)
            a[0] += np.where(viable_mask, np.nan_to_num(v), 0.0)
            a[1] += viable_mask.astype(float)
        for ph, (summ, nvi) in acc.items():
            with np.errstate(invalid="ignore"):
                per_lvl = np.where(nvi > 0, summ / nvi * 100.0, np.nan)
            mean_mut = np.nanmean(per_lvl)
            mean_n = np.mean(nvi)
            if mean_n >= 3:  # skip phyla with too few partners
                pw.writerow([strain, ph, f"{mean_mut:.3f}", f"{mean_n:.1f}"])
    print(f"[{sys_key}] phylum -> {ppath}")

if __name__ == "__main__":
    phylum = load_phylum()
    for k in (sys.argv[1:] or ["lac","akk"]):
        run(k, phylum)
