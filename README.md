# fast-mic-manuscript

**Reproducibility repository for the fast-mic manuscript** — exhaustive,
pair-resolved metabolic interaction typing between probiotics (*Akkermansia*,
*Lactobacillus*-group) and the human gut microbiome (UHGG, 3,238 species
representatives) across a ten-step prebiotic gradient (L0–L9).

> Associated study: *"fast-mic: a scalable tool for exhaustive pairwise
> interaction typing of genome-scale metabolic models reveals that carbon
> quality shapes probiotic–microbiome cooperation across a prebiotic gradient."*
>
> The fast-mic engine (Rust source) lives in a separate repository:
> <https://github.com/Zhangjyfree/fast-mic.git>. **This repository ships the
> data, benchmark inputs, and scripts** needed to reproduce every figure.

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
│   ├── run_gradient.sh                    # how the raw per-pair results were produced (needs the engine repo)
│   └── benchmark/
│       ├── benchmark_cobra.py             # COBRApy comparison
│       └── run_thread_scaling.sh          # thread-scaling harness
├── results/
│   ├── figures_paper/                     # 18 final figures × 3 formats
│   ├── akk_vs_uhgg/  lac_vs_uhgg/         # raw per-pair output, L0–L9
│   ├── fig1/                              # timing/memory tables + benchmark/
│   │   └── benchmark/                     # thread_scaling/ + correctness/ (COBRApy benchmark raw data)
│   ├── fig2/ fig4/ fig5/ fig6/            # main-figure intermediates
│   ├── figS1/ figS2/ figS4/ figS5/        # supplementary intermediates
│   ├── figS7/ figS9/ figS10/ figS11/ figS12/
│   ├── litvalidation/                     # literature cross-check
│   └── tableS5/cff_biomass_deviation.tsv  # 1,000 models × 10 media, FBA vs post-CFF biomass
└── test/
    ├── akk/  lac/                         # genomes, GEMs, IQ-TREE trees
    └── UHGG/                              # split UHGG model archive + GTDB metadata
```

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
S4**; the medium CSVs themselves live in the engine repository under `media/`
(see §8).

---

## 4. Requirements

**R (≥ 4.2)** — CRAN: `tidyverse` (ggplot2, dplyr, tidyr, readr, tibble),
`patchwork`, `scales`, `ggrepel`, `ggalluvial`, `cowplot`, `viridisLite`;
Bioconductor (Fig 2 only): `ggtree`, `treeio`, `aplot`.

```r
install.packages(c("tidyverse","patchwork","scales","ggrepel","ggalluvial","cowplot"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ggtree","treeio","aplot"))   # Fig 2 only
```

PDF output uses `cairo_pdf` so that Unicode glyphs (μ, β, →, ×) survive; a
`pdf()` device without cairo silently drops them.

**Python (≥ 3.8)** — the extract/table scripts use the standard library only,
except `extract_fig6_enrich.py` (needs `numpy` for the bootstrap) and
`stats_strain_level.py` (needs `numpy` and `scipy` for the exact rank tests).

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

**fast-mic binary** — only needed to re-run simulations:

```bash
git clone https://github.com/Zhangjyfree/fast-mic.git
cd fast-mic && cargo build --release        # → target/release/fast-mic
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
Fig 2 210, Fig 3 160, Fig 4 205, Fig 5 200, Fig 6 225, S1 105, S2 195, S3 70,
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

| Directory | Figure / table | Command | Inputs |
|---|---|---|---|
| `results/fig1/benchmark/` | 1B–1E | `bash scripts/benchmark/run_thread_scaling.sh` (COBRApy comparison + thread scaling; needs `cobra` and the fast-mic binary) | UHGG model corpora, L0–L9 media |
| `results/tableS5/cff_biomass_deviation.tsv` | Table S5 | run from the fast-mic checkout: `cargo build --release --bin bench-cff-deviation` then `./target/release/bench-cff-deviation --media-list media/gradient_media_list.txt --model-list <model list> > cff_biomass_deviation.tsv`. The per-model TSV goes to stdout and the validation summary (max / mean \|deviation\|, worst case) to stderr; Table S5 aggregates the TSV by prebiotic level. Both list files hold absolute paths — rewrite them for your own checkout. | fast-mic binary; the same 1,000 UHGG models as `results/fig1/benchmark/correctness/model_list.txt`, L0–L9 media |
| `results/fig1/` | 1F | `python3 scripts/scaling_timing_smetana.py --sys {akk,lac} --sizes 100,500,1000 --reps 3 --solver cplex` | fast-mic binary + SMETANA env |
| `results/fig1/…_mem.tsv` | 1G | same script with `--out-suffix _mem` — records `peak_rss_mb` (fast-mic per child process via `os.wait4`; SMETANA via `getrusage(RUSAGE_SELF)`) | fast-mic binary + SMETANA env |
| `results/fig2/` | 2 | `python3 scripts/extract_fig2_traits.py` | `test/{akk,lac}/…_gapseq_wdm_xml`, trees, `{akk,lac}_vs_uhgg/` |
| `results/fig4/` | 4G | `fast-mic` on L5 / L5-glc / L6 media (§8) | fast-mic binary |
| `results/fig5/` | 5 | `python3 scripts/extract_fig5_crossfeed.py` | `{akk,lac}_vs_uhgg/*.full.tsv`, ModelSEED `compounds.tsv` (§8) |
| `results/fig6/` | 6 | `python3 scripts/extract_fig6_enrich.py` | `{akk,lac}_vs_uhgg/`, `test/UHGG/…metadata…gz` (needs `numpy`) |
| `results/figS2/` | S2B,C | `python3 scripts/run_smetana_batch.py <level> <sys> 300` | `{akk,lac}_vs_uhgg/`, `test/` models, SMETANA env |
| `results/figS7/` | S7 | `python3 scripts/make_crossfeed_table.py` **and** `python3 scripts/extract_interception.py` | `{akk,lac}_vs_uhgg/*.full.tsv`, `compounds.tsv`, UHGG models |
| `results/figS9/` | S9 | `fast-mic` at L5/L6, default vs `--fixed-ratio` (§8) | fast-mic binary |
| `results/figS4/` | S4, Table S9 | `python3 scripts/stats_strain_level.py` (bootstrap CIs, exact Mann-Whitney and Wilcoxon; B = 10,000, seed 0) | `results/{akk,lac}_vs_uhgg/L5_pectin.tsv`, `L6_resistant_starch.tsv` |
| `results/litvalidation/` | Table S11 columns G–K | `python3 scripts/build_lit_validation_table.py` | `results/fig5/fig5_prevalence_{akk,lac}.tsv` + curated evidence in the script |
| `results/figS1/`, `figS5/`, `figS10/`–`figS12/` | S1, S5, S10–S12 | focused fast-mic runs / analyses — **provided**; regeneration needs the engine binary | — |

```bash
# Intermediates with a generator in this repository
python3 scripts/extract_fig2_traits.py        # → results/fig2/
python3 scripts/extract_fig6_enrich.py        # → results/fig6/   (numpy)
python3 scripts/extract_fig5_crossfeed.py     # → results/fig5/   (needs compounds.tsv, §8)
python3 scripts/make_crossfeed_table.py       # → results/figS7/crossfeed_landscape_table.tsv
python3 scripts/extract_interception.py       # → results/figS7/interception.tsv
python3 scripts/stats_strain_level.py         # → results/figS4/  (statistics for Fig S4 + Table S9)
python3 scripts/build_lit_validation_table.py # → results/litvalidation/lit_validation.tsv
```

**Literature cross-check** (`results/litvalidation/lit_validation.tsv`): every
predicted cross-feeding currency is scored against the experimental literature
as `supported` / `partial` / `untested` / `contradicted`, with a separate column
recording whether the published donor taxon matches the one predicted here.
Of the 27 currencies exported in at least a tenth of mutualistic pairs, 23 are
documented and one (*Lactobacillus* inosine export) points the other way.

---

## 7. Reassembling the UHGG model set

`test/UHGG/final_gapseq_xml/` ships the 3,238 gapseq SBML models as a tar
archive split into 44 × 20 MB parts (GitHub's per-file limit). Scripts that read
individual models (`extract_interception.py`, `run_smetana_batch.py`,
`run_gradient.sh`) expect the **extracted** `.xml` files:

```bash
cd test/UHGG/final_gapseq_xml
cat final_gapseq_xml.tar.gz.part_* > final_gapseq_xml.tar.gz
tar -xzf final_gapseq_xml.tar.gz          # → *.xml  (~3,238 models)
```

The extracted `.xml` files are git-ignored, so re-extracting never dirties the
working tree.

---

## 8. Not included here

| Item | Where | Why |
|---|---|---|
| fast-mic Rust source + binary | engine repo | separate software release |
| `media/gradient_L*_gapseq.csv` (12 files: L0–L9 + the L5-glc and L6′ controls) | engine repo, `media/` | composition and uptake bounds are tabulated in Table S4 |
| ModelSEED `compounds.tsv` | engine repo | compound id → name map, only needed to regenerate `results/fig5/` |
| COBRApy per-model raw logs | `results/fig1/benchmark/thread_scaling/` | 138 log/TSV files are included as-is |

`scripts/run_gradient.sh` is kept for provenance: it documents exactly how
`results/{akk,lac}_vs_uhgg/` was produced, but it needs the engine binary and
the medium CSVs, so it cannot run inside this repository as-is.

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
| Fig 4 | `plot_fig4_mechanisms.R` | A competition intensity · B Δmutualism vs Δcompetition · C/D growth · E cross-fed metabolites · F gene-supported fraction · G controlled pectin→glucose contrast |
| Fig 5 | `extract_fig5_crossfeed.py` → `plot_fig5_crossfeed.R` | A Sankey of currencies · B condition-resolved prevalence · C metabolite-class composition |
| Fig 6 | `extract_fig6_enrich.py` → `plot_fig6_heterogeneity.R` | A/B strain heatmaps · C bootstrap CI · D capacity does *not* predict cooperation · E partner phylum |
| Fig S1 | `plot_figS1_l0_substrates.R` | L0 substrate basis of *Lactobacillus* viability |
| Fig S2 | `run_smetana_batch.py` → `plot_figS2_smetana_comparison.R` | A projected cost · B potential vs realized (AUROC ≈ 0.5) · C rank scatters |
| Fig S3 | `plot_figS3_interaction_composition.R` | full six-category interaction composition |
| Fig S4 | `stats_strain_level.py` → `plot_figS4_strain_bootstrap.R` | A generalist vs specialist at L5 · B paired L5→L6 decline; every statistic is read from `strain_level_tests.tsv`, none is hard-coded |
| Fig S5 | `plot_figS5_benefit_metric.R` | relative vs absolute benefit metric |
| Fig S6 | `plot_figS6_benefit_landscape.R` | joint β_A–β_B density |
| Fig S7 | `make_crossfeed_table.py` + `extract_interception.py` → `plot_figS7_crossfeed_landscape.R` | A currency landscape · B higher-order interception bound |
| Fig S8 | `plot_figS8_threshold_sensitivity.R` | A viability threshold · B benefit threshold |
| Fig S9 | `plot_figS9_objective_sensitivity.R` | lexicographic vs fixed-ratio objective |
| Fig S10 | `plot_figS10_uhgg_sensitivity.R` | genome-set sensitivity (quality, CPR) |
| Fig S11 | `plot_figS11_gpr_coverage.R` | GPR coverage vs cross-feeding gene support |
| Fig S12 | `plot_figS12_gapfill_knockout.R` | A overall mutualism · B fate of mutualistic pairs after gap-fill knockout |

---

## 11. Citation

If you use these data or scripts, please cite the fast-mic manuscript.

- Engine: <https://github.com/Zhangjyfree/fast-mic.git>
- Reproduction: <https://github.com/Zhangjyfree/fast-mic-manuscript.git>
