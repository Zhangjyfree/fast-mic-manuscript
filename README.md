# fast-mic-manuscript

**Reproducibility repository for the fast-mic manuscript** — exhaustive,
pair-resolved metabolic interaction typing between probiotics (*Akkermansia*,
*Lactobacillus*-group) and the human gut microbiome (UHGG, 3,238 species
representatives) across a ten-step prebiotic gradient (L0–L9).

> Associated study: *"Fast-mic: a scalable tool for exhaustive pairwise
> interaction typing of genome-scale metabolic models reveals that carbon
> quality shapes probiotic–microbiome cooperation across a prebiotic gradient."*
>
> The fast-mic engine (Rust source) lives in a separate repository:
> <https://github.com/Zhangjyfree/fast-mic.git>. **This repository ships the
> data, benchmark inputs, and scripts** needed to reproduce every figure; the
> media are part of the engine repository (`media/`).

**Paths.** Run every script from the root of this repository; every path in
the scripts, model lists and logs is relative to it. Scripts that need the
engine (binary or media) look for it in `../fast-mic` (clone the two
repositories side by side) or wherever `FASTMIC_ENGINE` points.

---

## 1. What is here

| | |
|---|---|
| **6 main figures + 12 supplementary figures** | `results/figures_paper/` — each as `pdf` / `png` / `tiff` |
| **Every plotting script** | `scripts/plot_fig*.R` — one script per figure, no hidden steps |
| **Every intermediate table** | `results/fig*/`, `results/figS*/`, `results/litvalidation/` |
| **Raw per-pair simulation output** | `results/{akk,lac}_vs_uhgg/` — 10 media × all pairs |
| **Benchmark inputs** | `results/fig1/benchmark/` — COBRApy correctness + thread scaling; generators in `scripts/benchmark/` |
| **Genomes, models, trees** | `test/` — the two probiotic panels and the UHGG model set |

Only the **submitted version** of each figure and script is kept here.
Superseded drafts live in `_superseded/` in the working tree and are not
mirrored into this repository.

---

## 2. Directory layout

```
fast-mic-manuscript/
├── scripts/
│   ├── figure_theme.R                     # shared figure conventions: Arial, BMC geometry, italic taxon names
│   ├── plot_fig1_benchmark.R              # Fig 1 benchmark, panels A–G
│   ├── extract_fig2_traits.py             # Fig 2 GEM traits + per-level monoculture growth
│   ├── plot_fig2_repertoire.R             # Fig 2 phylogeny + model repertoire + growth heatmap
│   ├── plot_fig3_gradient_overview.R      # Fig 3 interactions across the gradient
│   ├── plot_fig4_mechanisms.R             # Fig 4 glucose-crash mechanisms, panels A–G
│   ├── extract_fig5_crossfeed.py          # Fig 5 per-level cross-feeding prevalence
│   ├── plot_fig5_crossfeed.R              # Fig 5 Sankey + condition heatmaps + class composition
│   ├── extract_fig6_enrich.py             # Fig 6 bootstrap CIs + partner phyla
│   ├── plot_fig6_heterogeneity.R          # Fig 6 strain heterogeneity, panels A–E
│   ├── plot_figS1_l0_substrates.R … plot_figS12_gapfill_knockout.R   # Supplementary Figs S1–S12
│   ├── scaling_timing_smetana.py          # fast-mic vs SMETANA timing + peak memory (Fig 1F/1G)
│   ├── run_smetana_batch.py               # SMETANA MIP/MRO batch runs (Fig S2B/C)
│   ├── extract_interception.py            # higher-order interception bound (Fig S7B)
│   ├── make_crossfeed_table.py            # pooled cross-feeding currency table (Fig S7A)
│   ├── stats_strain_level.py              # strain-level statistics: bootstrap CIs, Mann-Whitney, Wilcoxon (Fig S4 + Table S9)
│   ├── build_lit_validation_table.py      # predicted currencies vs published experimental evidence
│   ├── figS1_l0_uptake.py                 # L0 uptake by substrate class (Fig S1, Table S3; engine)
│   ├── figS5_benefit_metric.py            # relative vs absolute benefit (Fig S5)
│   ├── figS9_objective_sensitivity.py     # lexicographic vs fixed-ratio objective (Fig S9, Table S12; engine)
│   ├── figS10_uhgg_sensitivity.py         # UHGG genome-set sensitivity (Fig S10)
│   ├── figS11_gpr_coverage.py             # gene support of reactions vs cross-feeding flux (Fig S11)
│   ├── gapfill_knockout.py                # gap-fill-candidate knockout (Fig S12; engine)
│   ├── tableS5_summary.py                 # loop-removal validation by level (Table S5)
│   ├── run_fig4_contrast.sh               # controlled pectin/glucose contrast (Fig 4G, Table S10; engine)
│   ├── akk_sugar_block.py                 # Akkermansia free-sugar uptake blocked (Table S10, Discussion; engine)
│   ├── run_gradient.sh                    # the L0–L9 screen behind results/{akk,lac}_vs_uhgg/ (engine)
│   └── benchmark/
│       ├── benchmark_cobra.py             # COBRApy comparison
│       ├── accuracy_check.py              # fast-mic vs COBRApy agreement statistics
│       └── run_thread_scaling.sh          # thread-scaling + correctness harness (engine)
├── results/
│   ├── figures_paper/                     # 18 final figures × 3 formats
│   ├── akk_vs_uhgg/  lac_vs_uhgg/         # raw per-pair output, L0–L9
│   ├── fig1/                              # timing/memory tables + benchmark/
│   │   └── benchmark/                     # thread_scaling/ + correctness/ (COBRApy benchmark raw data)
│   ├── fig2/ fig4/ fig5/ fig6/            # main-figure intermediates
│   ├── figS1/ figS2/ figS4/ figS5/        # supplementary intermediates
│   ├── figS7/ figS9/ figS10/ figS11/ figS12/
│   ├── litvalidation/                     # literature cross-check
│   ├── tableS5/                           # 1,000 models × 10 media, FBA vs post-CFF biomass (+ by-level summary)
│   └── tableS10/                          # Akkermansia free-sugar uptake control (Table S10, last three rows)
└── test/
    ├── akk/  lac/                         # genomes, GEMs, IQ-TREE trees
    └── UHGG/                              # split UHGG model archive + GTDB metadata
```

**Which probiotic models were used.** Every analysis in the paper uses
`test/akk/akk_gapseq_xml/` (*Akkermansia*, gap-filled on the *Akkermansia*
minimal medium) and `test/lac/lac_genomes_faa_gapseq_wdm_xml/`
(*Lactobacillus*-group, gap-filled on Western diet + mucin), as stated in the
Methods.

---

## 3. Prebiotic gradient (L0–L9)

Ten strictly nested, cumulative media. Each level adds the *in vivo* hydrolysis
products of one more prebiotic on top of all previous levels. **L6 is the first
level to release free glucose** (the "glucose crash").

| Level | File                  | Added carbon                                              |
|-------|-----------------------|-----------------------------------------------------------|
| L0    | `L0_base`             | housekeeping nutrients + mucin glycan monomers            |
| L1    | `L1_inulin`           | inulin (D-fructose, sucrose)                              |
| L2    | `L2_fos`              | FOS (inulobiose)                                          |
| L3    | `L3_gos`              | GOS (lactose, lactulose, D-galactose)                     |
| L4    | `L4_xos`              | XOS (D-xylose, L-arabinose)                               |
| L5    | `L5_pectin`           | pectin (D-galacturonate, L-rhamnose) — peak               |
| L6    | `L6_resistant_starch` | resistant starch (D-glucose, maltose) — **free glucose**  |
| L7    | `L7_bglucan`          | β-glucan (cellobiose)                                     |
| L8    | `L8_hmo`              | human milk oligosaccharides (lacto-N-biose)               |
| L9    | `L9_mos`              | MOS (D-mannose, mannobiose)                               |

Medium composition and uptake bounds are tabulated in **Supplementary Table
S4**; the medium CSVs are in the engine repository under `media/` (the ten
levels are listed, in order, in `media/gradient_media_list.txt`). The same
directory holds the L5-glc control of Figure 4G (`gradient_L5glc_gapseq.csv`:
the L1–L4 background with pectin replaced by free glucose) and the two
gap-filling media (`akk_minimal_gapseq.csv` for *Akkermansia*,
`western_diet_mucin_gapseq.csv` for the other models).

---

## 4. Requirements

**R (≥ 4.2)** — CRAN: `tidyverse` (ggplot2, dplyr, tidyr, readr, tibble),
`patchwork`, `scales`, `ggrepel`, `ggalluvial`, `cowplot`, `viridisLite`, `ggtext`;
Bioconductor (Fig 2 only): `ggtree`, `treeio`, `aplot`.

```r
install.packages(c("tidyverse","patchwork","scales","ggrepel","ggalluvial","cowplot","ggtext"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ggtree","treeio","aplot"))   # Fig 2 only
```

PDF output uses `cairo_pdf` so that Unicode glyphs (μ, β, →, ×) survive; a
`pdf()` device without cairo silently drops them.

**Python (≥ 3.8)** — the extract/table scripts use the standard library only,
except `extract_fig6_enrich.py` (needs `numpy` for the bootstrap),
`stats_strain_level.py` (needs `numpy` and `scipy` for the exact rank tests) and
`gapfill_knockout.py` and `scripts/benchmark/` (need `cobra`).

**SMETANA comparison only** (Fig 1F, 1G, S2) — a conda env with `smetana`,
`reframed` and a commercial MILP solver:

```bash
conda create -n smetana python=3.10 && conda activate smetana
pip install smetana reframed        # + CPLEX 22.2 (or Gurobi) Python bindings
```

`reframed` has **no native HiGHS interface** (its only HiGHS route is a PuLP
wrapper that the default solver order `gurobi > cplex > scip` never selects), so
SMETANA is run on CPLEX while fast-mic keeps its free, natively linked HiGHS.
The asymmetry favours SMETANA and is stated as such in the manuscript.

**fast-mic engine** — only needed to re-run simulations. Clone it next to this
repository (or set `FASTMIC_ENGINE` to its location); it also provides the
media (`media/gradient_*_gapseq.csv`) and the ModelSEED `media/compounds.tsv`
used by `extract_fig5_crossfeed.py` and `make_crossfeed_table.py`:

```bash
git clone https://github.com/Zhangjyfree/fast-mic.git ../fast-mic
(cd ../fast-mic && cargo build --release)   # → ../fast-mic/target/release/fast-mic
```

---

## 5. Reproduce the figures

Every plotting script sources `scripts/figure_theme.R`, which holds the
conventions the journal expects: Arial throughout, BMC page geometry
(170 mm full-page / 85 mm half-page, 225 mm maximum height, 300 dpi, lines
above 0.25 pt) and italic genus/species names. The helpers are `sp_md()` for
scale and label text, `sp_labeller()` for facet strips and `sp_plotmath()` for
geoms that cannot render markdown (ggtree tip labels). They change only how a
label is drawn — never the underlying data — and the text element that shows
them is an `element_markdown()`, so run the scripts from the repository root
so that `scripts/figure_theme.R` resolves.

Every figure is drawn at its final printed size: 170 mm wide (BMC full page),
at most 225 mm tall, 300 dpi, Arial throughout. Heights by figure: Fig 1 225,
Fig 2 210, Fig 3 160, Fig 4 225, Fig 5 200, Fig 6 225, S1 105, S2 195, S3 70,
S4 82, S5 88, S6 150, S7 190, S8 80, S9 95, S10 84, S11 80, S12 85 mm.

```bash
# Main figures
# Fig 1 benchmark data sits in results/fig1/benchmark/ — no archive to unpack
Rscript scripts/plot_fig1_benchmark.R      # default input paths point there
Rscript scripts/plot_fig2_repertoire.R    results/figures_paper/fig2_repertoire
Rscript scripts/plot_fig3_gradient_overview.R
Rscript scripts/plot_fig4_mechanisms.R
Rscript scripts/plot_fig5_crossfeed.R     results/figures_paper/fig5_crossfeed
Rscript scripts/plot_fig6_heterogeneity.R results/figures_paper/fig6_heterogeneity

# Supplementary figures S1–S12 (numbered in order of first citation)
Rscript scripts/plot_figS1_l0_substrates.R           --outdir results/figures_paper
Rscript scripts/plot_figS2_smetana_comparison.R      results/figures_paper/figS2_smetana_comparison
Rscript scripts/plot_figS3_interaction_composition.R
Rscript scripts/plot_figS4_strain_bootstrap.R        --outdir results/figures_paper
Rscript scripts/plot_figS5_benefit_metric.R          --outdir results/figures_paper
Rscript scripts/plot_figS6_benefit_landscape.R
Rscript scripts/plot_figS7_crossfeed_landscape.R     # needs the two tables in §6
Rscript scripts/plot_figS8_threshold_sensitivity.R
Rscript scripts/plot_figS9_objective_sensitivity.R
Rscript scripts/plot_figS10_uhgg_sensitivity.R       --outdir results/figures_paper
Rscript scripts/plot_figS11_gpr_coverage.R           --outdir results/figures_paper
Rscript scripts/plot_figS12_gapfill_knockout.R       --outdir results/figures_paper
```

Every intermediate table these scripts read is already in `results/`, so the
figures reproduce without re-running any simulation.

---

## 6. Intermediate files — how each is generated

All commands run from the repository root. "Engine" means the fast-mic binary
and media from `../fast-mic` (or `$FASTMIC_ENGINE`); commands that read
individual UHGG models need the extracted model set (§7).

| Directory | Figure / table | Command | Inputs |
|---|---|---|---|
| `results/akk_vs_uhgg/`, `results/lac_vs_uhgg/` | 3, 5, 6, S3–S8, S10, Tables S8, S9, S11 | `bash scripts/run_gradient.sh --group1 test/akk/akk_gapseq_xml --group2 test/UHGG/final_gapseq_xml --threads 12 --full-tsv --out results/akk_vs_uhgg` (and the same with `test/lac/lac_genomes_faa_gapseq_wdm_xml` → `results/lac_vs_uhgg`) | engine |
| `results/fig1/benchmark/` | 1B–1E | `N_LIST=100,500,1000 CORRECTNESS_N=1000 bash scripts/benchmark/run_thread_scaling.sh` (COBRApy comparison + thread scaling; builds `bench-single-fba` in the engine checkout) | engine, `cobra`, UHGG models |
| `results/tableS5/` | Table S5 | `(cd ../fast-mic && cargo build --release --bin bench-cff-deviation)`, then `../fast-mic/target/release/bench-cff-deviation --media-list results/fig1/benchmark/media_list.txt --model-list results/fig1/benchmark/correctness/model_list.txt > results/tableS5/cff_biomass_deviation.tsv 2> results/tableS5/cff_biomass_deviation_summary.txt`, then `python3 scripts/tableS5_summary.py` (by-level table) | engine, the 1,000 correctness models |
| `results/fig1/` | 1F | `python3 scripts/scaling_timing_smetana.py --sys {akk,lac} --sizes 100,500,1000 --reps 3 --solver cplex --smetana-cap-pairs 2000` (the cap must exceed 1,002, the largest *Akkermansia* size: 6 strains × 167 UHGG genomes) | engine + SMETANA env |
| `results/fig1/…_mem.tsv` | 1G | same script with `--reps 1 --out-suffix _mem` — records `peak_rss_mb` (fast-mic per child process via `os.wait4`; SMETANA via `getrusage(RUSAGE_SELF)`) | engine + SMETANA env |
| `results/fig2/` | 2, 6D | `python3 scripts/extract_fig2_traits.py` | `test/akk/akk_gapseq_xml`, `test/lac/lac_genomes_faa_gapseq_wdm_xml`, trees, `{akk,lac}_vs_uhgg/` |
| `results/fig4/` | 4G, Table S10 | `bash scripts/run_fig4_contrast.sh 12` | engine (media incl. L5-glc) |
| `results/tableS10/` | Table S10 (last three rows), Discussion | `python3 scripts/akk_sugar_block.py --threads 12` — blocks the *Akkermansia* glucose/maltose transporters (no gene association in any strain), checks that its monoculture growth is then identical on L5, L5-glc and L6, and re-screens the three media | engine |
| `results/fig5/` | 5 | `python3 scripts/extract_fig5_crossfeed.py` | `{akk,lac}_vs_uhgg/*.full.tsv[.gz]`, engine `media/compounds.tsv` |
| `results/fig6/` | 6 | `python3 scripts/extract_fig6_enrich.py` | `{akk,lac}_vs_uhgg/`, `test/UHGG/…metadata…gz` (needs `numpy`) |
| `results/figS1/` | S1, Table S3 | `python3 scripts/figS1_l0_uptake.py` | engine, `test/lac/…` |
| `results/figS2/` | S2B,C | `python3 scripts/run_smetana_batch.py <L5_pectin\|L6_resistant_starch> <akk\|lac> 300` | `{akk,lac}_vs_uhgg/`, `test/` models, SMETANA env |
| `results/figS4/` | S4, Table S9 | `python3 scripts/stats_strain_level.py` (bootstrap CIs, exact Mann-Whitney and Wilcoxon; B = 10,000, seed 0) | `results/{akk,lac}_vs_uhgg/L5_pectin.tsv`, `L6_resistant_starch.tsv` |
| `results/figS5/` | S5 | `python3 scripts/figS5_benefit_metric.py` | `{akk,lac}_vs_uhgg/` |
| `results/figS7/` | S7 | `python3 scripts/make_crossfeed_table.py` **and** `python3 scripts/extract_interception.py` | `{akk,lac}_vs_uhgg/*.full.tsv[.gz]`, engine `media/compounds.tsv`, UHGG models |
| `results/figS9/` | S9, Table S12 | `python3 scripts/figS9_objective_sensitivity.py --threads 12` (*L. gasseri* × the 3,218 bacterial UHGG models; the 20 archaea are excluded via the UHGG metadata) | engine |
| `results/figS10/` | S10 | `python3 scripts/figS10_uhgg_sensitivity.py` | `{akk,lac}_vs_uhgg/`, `test/UHGG/…metadata…gz` |
| `results/figS11/` | S11 | `python3 scripts/figS11_gpr_coverage.py` | probiotic models, `{akk,lac}_vs_uhgg/` |
| `results/figS12/R15_summary.tsv` | S12 | `python3 scripts/gapfill_knockout.py --system Akkermansia --models test/akk/akk_gapseq_xml --partners test/UHGG/final_gapseq_xml --control results/akk_vs_uhgg/L5_pectin.tsv --threads 12 --append`, then the same with `--system Lactobacillus --models test/lac/lac_genomes_faa_gapseq_wdm_xml --control results/lac_vs_uhgg/L5_pectin.tsv` | engine, `cobra` |
| `results/litvalidation/` | Table S11 columns G–K | `python3 scripts/build_lit_validation_table.py` | `results/fig5/fig5_prevalence_{akk,lac}.tsv` + curated evidence in the script |

```bash
# Intermediates computed from the stored per-pair tables (no engine needed)
python3 scripts/extract_fig2_traits.py        # → results/fig2/
python3 scripts/extract_fig6_enrich.py        # → results/fig6/   (numpy)
python3 scripts/extract_fig5_crossfeed.py     # → results/fig5/   (engine media/compounds.tsv)
python3 scripts/make_crossfeed_table.py       # → results/figS7/crossfeed_landscape_table.tsv
python3 scripts/extract_interception.py       # → results/figS7/interception.tsv (UHGG models)
python3 scripts/stats_strain_level.py         # → results/figS4/  (statistics for Fig S4 + Table S9)
python3 scripts/build_lit_validation_table.py # → results/litvalidation/lit_validation.tsv
python3 scripts/figS5_benefit_metric.py       # → results/figS5/
python3 scripts/figS10_uhgg_sensitivity.py    # → results/figS10/
python3 scripts/figS11_gpr_coverage.py        # → results/figS11/
python3 scripts/tableS5_summary.py            # → results/tableS5/tableS5_by_level.tsv
```

Every script above reproduces the stored file byte for byte. Re-running an
engine command reproduces the stored interaction types exactly; floating-point
columns can differ in the last printed digits (≤ 1×10⁻⁵) because multi-threaded
LP solves are not bit-reproducible. `run_gradient.sh` takes several hours per
system at 12 threads; `run_fig4_contrast.sh` about 40 min.

**Literature cross-check** (`results/litvalidation/lit_validation.tsv`): every
predicted cross-feeding currency is scored against the experimental literature
as `supported` / `partial` / `untested` / `contradicted`, with a separate column
recording whether the published donor taxon matches the one predicted here.

The file holds all 33 curated currencies, while the manuscript quotes the 27
that pass a threshold: **mean** export prevalence across L0–L9 of at least 10%
("exported in at least a tenth of mutualistic pairs"). The cut is on the mean,
not on the per-level maximum — *Lactobacillus* methanethiol and 1,2-propanediol
reach 10% at their best level but not on average, so a max-based cut would give
29 currencies instead of 27. Of those 27, 23 are documented (17 directly, 6
partly) and one (*Lactobacillus* inosine export) points the other way. The
script prints both counts, so the figures quoted in the text can be checked
without re-deriving the filter:

```
all curated currencies: contradicted=1, partial=6, supported=22, untested=4
mean export prevalence >= 10% (the cut quoted in the manuscript): n=27, contradicted=1, partial=6, supported=17, untested=3
  -> 23 of 27 documented (17 directly, 6 partly); 3 untested, 1 contradicted
```

---

## 7. Reassembling the UHGG model set

`test/UHGG/final_gapseq_xml/` ships the 3,238 gapseq SBML models as a tar
archive split into 44 × 20 MB parts (GitHub's per-file limit). Scripts that read
individual models (`run_gradient.sh`, `run_fig4_contrast.sh`, `gapfill_knockout.py`,
`figS1_l0_uptake.py`, `figS9_objective_sensitivity.py`, `extract_interception.py`,
`run_smetana_batch.py`, `scaling_timing_smetana.py`, the `benchmark/` harness and
Table S5) expect the **extracted** `.xml` files (about 17 GB):

```bash
# from the repository root; the archive holds a top-level final_gapseq_xml/ directory
cat test/UHGG/final_gapseq_xml/final_gapseq_xml.tar.gz.part_* | tar -xzf - -C test/UHGG
ls test/UHGG/final_gapseq_xml/*.xml | wc -l   # → 3238
```

The extracted `.xml` files are git-ignored, so re-extracting never dirties the
working tree.

---

## 8. Not included here

| Item | Where | Why |
|---|---|---|
| fast-mic Rust source + binary, `bench-single-fba`, `bench-cff-deviation` | engine repo | separate software release |
| Media (`gradient_L*_gapseq.csv` for L0–L9 and L5-glc, gap-filling media) | engine repo, `media/` | shared with the engine; composition in Table S4 |
| ModelSEED `compounds.tsv` | engine repo, `media/compounds.tsv` | compound id → name map (12 MB), read by `extract_fig5_crossfeed.py` and `make_crossfeed_table.py` |
| GTDB-Tk de novo alignments (input to IQ-TREE) | not included | the trees, model selection and support values are in `test/{akk,lac}/iqtree_out/`; each `*.log` records the full IQ-TREE 3.1.2 command (GTDB-Tk v2.7 `de_novo_wf`, R232) |

The COBRApy/fast-mic benchmark logs and model lists in `results/fig1/benchmark/`
and the IQ-TREE logs are the files written by the original runs; only the
absolute path prefixes they recorded were rewritten to the repository-relative
paths used here (`test/UHGG/final_gapseq_xml/`, `../fast-mic/media/`, `results/fig1/benchmark/`,
`iqtree3`).

---

## 9. Data dictionary

Each level produces a summary `L*.tsv` (one row per species pair). With
`--full-tsv`, a richer `L*.full.tsv` adds per-metabolite cross-feeding columns.
(`lac_vs_uhgg/*.full.tsv` are gzip-compressed and `akk_vs_uhgg/*.full.tsv` are
not; the readers handle both.)

| Column                        | Meaning                                                  |
|-------------------------------|----------------------------------------------------------|
| `species_a` / `species_b`     | the two members of the pair                              |
| `growth_a_alone` / `_b_alone` | monoculture growth rate (h⁻¹)                            |
| `growth_a_co` / `_b_co`       | co-culture growth rate (h⁻¹)                             |
| `benefit_a` / `benefit_b`     | relative benefit (co-culture vs monoculture)             |
| `interaction_type`            | mutualism, competition, commensalism, parasitism, …      |
| `competition_intensity`       | shared-uptake overlap                                    |
| `n_exchanged_metabolites`     | number of cross-fed metabolites                          |
| `gene_supported_fraction`     | gene-supported fraction of cross-feeding flux            |
| `a_to_b_metabolites` ⁺        | metabolites, fluxes and donor/receiver genes exchanged a→b |
| `competed_metabolites` ⁺      | metabolites competed for                                 |

⁺ `*.full.tsv` only.

A viable pair = both members grow in monoculture (> 1×10⁻⁴ h⁻¹); mutualism % is
computed over viable pairs.

The timing tables (`results/fig1/scaling_timing_*.tsv`) carry `tool`, `n_pairs`,
`n_prob`, `n_uhgg`, `threads`, `solver`, `rep`, `wall_s`, `sec_per_pair` and —
in the `_mem` runs — `peak_rss_mb`. Each size is a complete probiotic × UHGG
cross-product, so the realized pair counts are 100/500/1,000 for the ten-strain
*Lactobacillus* panel and 102/498/1,002 for the six-strain *Akkermansia* panel.

---

## 10. Figure index

| Figure | Script | Panels / topic |
|---|---|---|
| Fig 1 | `plot_fig1_benchmark.R` (+ `scaling_timing_smetana.py`) | A pipeline · B accuracy vs COBRApy (r = 1.000) · C runtime · D parallel speedup · E memory vs COBRApy · F interaction-typing scalability vs SMETANA · G interaction-typing memory vs SMETANA |
| Fig 2 | `extract_fig2_traits.py` → `plot_fig2_repertoire.R` | A *Akkermansia* · B *Lactobacillus*: phylogeny + model repertoire + monoculture growth (shared colour scale) |
| Fig 3 | `plot_fig3_gradient_overview.R` | A viable pairs · B mutualism · C competition · D mean net benefit |
| Fig 4 | `run_fig4_contrast.sh` (G; `akk_sugar_block.py` for the Table S10 control) → `plot_fig4_mechanisms.R` | A competition intensity · B Δmutualism vs Δcompetition · C/D growth · E cross-fed metabolites · F gene-supported fraction · G controlled pectin→glucose contrast |
| Fig 5 | `extract_fig5_crossfeed.py` → `plot_fig5_crossfeed.R` | A Sankey of currencies · B condition-resolved prevalence · C metabolite-class composition |
| Fig 6 | `extract_fig6_enrich.py` → `plot_fig6_heterogeneity.R` | A/B strain heatmaps · C bootstrap CI · D capacity does *not* predict cooperation · E partner phylum |
| Fig S1 | `figS1_l0_uptake.py` → `plot_figS1_l0_substrates.R` | L0 substrate basis of *Lactobacillus* viability |
| Fig S2 | `run_smetana_batch.py` → `plot_figS2_smetana_comparison.R` | A projected cost · B potential vs realized (AUROC 0.53–0.56) · C rank scatters |
| Fig S3 | `plot_figS3_interaction_composition.R` | full six-category interaction composition |
| Fig S4 | `stats_strain_level.py` → `plot_figS4_strain_bootstrap.R` | A generalist vs specialist at L5 · B paired L5→L6 decline; every statistic is read from `strain_level_tests.tsv`, none is hard-coded |
| Fig S5 | `figS5_benefit_metric.py` → `plot_figS5_benefit_metric.R` | relative vs absolute benefit metric |
| Fig S6 | `plot_figS6_benefit_landscape.R` | joint β_A–β_B density |
| Fig S7 | `make_crossfeed_table.py` + `extract_interception.py` → `plot_figS7_crossfeed_landscape.R` | A currency landscape · B higher-order interception bound |
| Fig S8 | `plot_figS8_threshold_sensitivity.R` | A viability threshold · B benefit threshold |
| Fig S9 | `figS9_objective_sensitivity.py` → `plot_figS9_objective_sensitivity.R` | lexicographic vs fixed-ratio objective |
| Fig S10 | `figS10_uhgg_sensitivity.py` → `plot_figS10_uhgg_sensitivity.R` | genome-set sensitivity (quality, CPR) |
| Fig S11 | `figS11_gpr_coverage.py` → `plot_figS11_gpr_coverage.R` | GPR coverage vs cross-feeding gene support |
| Fig S12 | `gapfill_knockout.py` → `plot_figS12_gapfill_knockout.R` | A overall mutualism · B fate of mutualistic pairs after gap-fill knockout |

---

## 11. Citation

If you use these data or scripts, please cite the fast-mic manuscript.

- Engine: <https://github.com/Zhangjyfree/fast-mic.git>
- Reproduction: <https://github.com/Zhangjyfree/fast-mic-manuscript.git>
