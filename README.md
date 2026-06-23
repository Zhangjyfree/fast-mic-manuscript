# fast-mic-manuscript

**Reproducibility repository for the fast-mic manuscript** — pairwise metabolic
interactions between probiotics (*Akkermansia*, *Lactobacillus*-group) and the
human gut microbiome (UHGG) across a 10-step prebiotic gradient (L0–L9).

**fast-mic 论文的复现仓库** —— 益生菌（*Akkermansia*、*Lactobacillus* 类群）
与人体肠道菌群（UHGG）之间的成对代谢互作，跨越 10 级益生元梯度（L0–L9）。

> Associated study / 关联研究: *"Genome-scale modelling indicates that carbon
> quality governs how generalist and specialist probiotics cooperate with the
> resident microbiome across a prebiotic gradient."*
> The fast-mic engine (Rust source) lives in a separate repository:
> <https://github.com/Zhangjyfree/fast-mic.git>. **This repository ships only the
> data + scripts** needed to reproduce the figures and tables. / fast-mic 引擎
> （Rust 源码）位于独立仓库；**本仓库只提供复现图表所需的数据 + 脚本**。

---

## 1. Overview / 概览

**EN** — This repo bundles the raw fast-mic simulation outputs, the genome /
model inputs, and every plotting script needed to regenerate the manuscript
figures: **main Figures 1–6** and **Supplementary Figures S1–S8**. The
provided outputs in `results/` let you redraw every figure directly; you only
need the compiled `fast-mic` binary (from the engine repo) to *re-run* the
underlying simulations.

**中文** —— 本仓库收录 fast-mic 的原始模拟输出、基因组/模型输入，以及复现**正文图
1–6** 和**补充图 S1–S8** 所需的全部绘图脚本。`results/` 中已附带输出，可直接重绘所
有图；只有在**重跑模拟**时才需要从引擎仓库编译得到 `fast-mic` 可执行文件。

---

## 2. Directory layout / 目录结构

```
fast-mic-manuscript/
├── README.md
├── scripts/                            # 全部脚本 / all scripts
│   ├── run_gradient.sh                 # 跑 fast-mic 梯度模拟 / run the gradient simulation
│   ├── make_crossfeed_table.py         # 汇总交叉喂养表 / pool the cross-feeding table
│   ├── plot_fig1_benchmark.R           # 图1 性能基准 / Fig 1 benchmark
│   ├── plot_fig2_phylo_trees.R         # 图2 系统发育树 / Fig 2 phylogeny
│   ├── plot_fig3_gradient_overview.R   # 图3 梯度总览 / Fig 3 gradient overview
│   ├── plot_fig4_mechanisms.R          # 图4 机制 / Fig 4 mechanisms
│   ├── plot_fig5_crossfeed_sankey.R    # 图5 交叉喂养桑基图 / Fig 5 Sankey
│   ├── plot_fig6_strain_heatmap.R      # 图6 菌株热图 / Fig 6 strain heatmap
│   └── plot_figS1…S8_*.R               # 补充图 S1–S8 / Supplementary figures (见第 8 节)
├── results/                            # fast-mic 输出 + 成图 / outputs + figures
│   ├── akk_vs_uhgg/                    # Akkermansia × UHGG, L0–L9 (.tsv + .full.tsv)
│   ├── lac_vs_uhgg/                    # Lactobacillus × UHGG, L0–L9 (.tsv + .full.tsv.gz)
│   ├── l0_substrates/                  # L0 单培养底物摄取 / L0 uptake by substrate class (Fig S1)
│   ├── glucose_isolation/             # 受控葡萄糖对照 / controlled glucose contrast (Fig S7)
│   ├── objective_sensitivity/          # fixed-ratio vs lexicographic (Fig S6)
│   ├── gpr_analysis/                   # 交叉喂养基因支持 vs GPR 覆盖 / GPR coverage (Fig S8)
│   └── figures_paper/                  # 最终 PDF/PNG/TIFF + 汇总表 / final figures + tables
└── test/                               # 输入基因组与模型 / input genomes & models
    ├── akk/                            # Akkermansia 基因组 / gapseq 模型 / 树
    ├── lac/                            # Lactobacillus 基因组 / gapseq 模型 / 树
    └── UHGG/                           # 肠道群落模型 + 元数据 / community models + metadata
```

**Two systems / 两个系统**

| Key 键        | Probiotic 益生菌                     | Community 群落 |
|---------------|--------------------------------------|----------------|
| `akk_vs_uhgg` | *Akkermansia* (6 strains, 3 species) | UHGG (gut)     |
| `lac_vs_uhgg` | *Lactobacillus*-group (10 strains)   | UHGG (gut)     |

> `akk_vs_uhgg` = specialist × gut; `lac_vs_uhgg` = generalist × gut. The
> specialist/generalist contrast runs through every figure. / `akk` 为特化菌、
> `lac` 为泛能菌；特化 vs 泛能的对照贯穿全部图表。

---

## 3. Prebiotic gradient / 益生元梯度 (L0–L9)

Ten strictly nested, cumulative media. Each level adds the *in vivo* hydrolysis
products of one more prebiotic on top of all previous levels. **L6 is the first
level to release free glucose** ("glucose crash"). / 十级严格嵌套、累加式培养基；
每一级在前面所有级别基础上加入一种益生元的体内水解产物；**L6 是首个释放游离
葡萄糖的级别**（“葡萄糖崩溃”）。

| Level | File 文件             | Added carbon 新增碳源 / 益生元                              |
|-------|-----------------------|-----------------------------------------------------------|
| L0    | `L0_base`             | 仅基础营养 + 黏蛋白聚糖（单体）/ housekeeping + mucin glycan monomers |
| L1    | `L1_inulin`           | 菊粉 inulin (D-fructose, sucrose)                         |
| L2    | `L2_fos`              | 低聚果糖 FOS (inulobiose)                                 |
| L3    | `L3_gos`              | 低聚半乳糖 GOS (lactose, lactulose, D-galactose)          |
| L4    | `L4_xos`              | 低聚木糖 XOS (D-xylose, L-arabinose)                      |
| L5    | `L5_pectin`           | 果胶 pectin (D-galacturonate, L-rhamnose) — 峰值 peak     |
| L6    | `L6_resistant_starch` | 抗性淀粉 (D-glucose, maltose) — **游离葡萄糖 free glucose**|
| L7    | `L7_bglucan`          | β-葡聚糖 (cellobiose)                                     |
| L8    | `L8_hmo`              | 母乳低聚糖 HMO (lacto-N-biose)                            |
| L9    | `L9_mos`              | 甘露寡糖 MOS (D-mannose, mannobiose)                      |

---

## 4. Requirements / 依赖

**fast-mic binary** — build from the engine repo (only needed to re-run
simulations): / fast-mic 可执行文件（仅重跑模拟时需要）：

```bash
git clone https://github.com/Zhangjyfree/fast-mic.git
cd fast-mic && cargo build --release    # → target/release/fast-mic
```

**R (≥ 4.2)** — CRAN: `tidyverse` (ggplot2, dplyr, tidyr, readr, tibble),
`patchwork`, `scales`, `ggrepel`, `ggalluvial`; Bioconductor (Fig 2 only):
`ggtree`, `treeio`.

```r
install.packages(c("tidyverse","patchwork","scales","ggrepel","ggalluvial"))
# Bioconductor (图2系统发育树 / Fig 2 phylogeny):
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ggtree","treeio"))
```

**Python (≥ 3.8)** — standard library only (`csv`, `os`, `collections`); no
third-party packages. / 仅用标准库，无需额外安装。

---

## 5. Reproduce the figures / 复现图表

**EN** — The simulation outputs in `results/` are already provided, so every
figure can be regenerated directly. Each R script defaults to reading from
`results/` and writing PDF + PNG + TIFF into `results/figures_paper/`; the input
and output paths are documented in the header comment of each script.

**中文** —— `results/` 中已附带模拟输出，可直接重绘所有图。每个 R 脚本默认从
`results/` 读取，并把 PDF + PNG + TIFF 输出到 `results/figures_paper/`；各脚本头部
注释标明了输入/输出路径。

```bash
# 正文图 / Main figures
Rscript scripts/plot_fig1_benchmark.R \
    benchmark/benchmark_results/thread_scaling/thread_scaling_results.tsv \
    benchmark/benchmark_results/correctness/scatter.tsv \
    benchmark/benchmark_results/correctness/stats.tsv     # ↑ benchmark data lives in the engine repo
Rscript scripts/plot_fig2_phylo_trees.R
Rscript scripts/plot_fig3_gradient_overview.R
Rscript scripts/plot_fig4_mechanisms.R
Rscript scripts/plot_fig5_crossfeed_sankey.R             # 需先生成下方汇总表 / needs the table below
Rscript scripts/plot_fig6_strain_heatmap.R

# 补充图 S1–S8 / Supplementary figures (citation order)
Rscript scripts/plot_figS1_l0_substrates.R               # in: results/l0_substrates/
Rscript scripts/plot_figS2_interaction_composition.R     # in: results/{akk,lac}_vs_uhgg/
Rscript scripts/plot_figS3_benefit_landscape.R
Rscript scripts/plot_figS4_crossfeed_landscape.R         # in: crossfeed_landscape_table.tsv (below)
Rscript scripts/plot_figS5_threshold_sensitivity.R       # in: results/{akk,lac}_vs_uhgg/
Rscript scripts/plot_figS6_objective_sensitivity.R       # in: results/objective_sensitivity/
Rscript scripts/plot_figS7_glucose_isolation.R           # in: results/glucose_isolation/
Rscript scripts/plot_figS8_gpr_coverage.R                # in: results/gpr_analysis/
```

**Cross-feeding summary table / 交叉喂养汇总表** (used by Fig 5 & Fig S4):

```bash
python3 scripts/make_crossfeed_table.py
# → results/figures_paper/crossfeed_landscape_table.tsv
```

It pools every mutualistic pair across L0–L9 and counts, for each exchanged
metabolite, the fraction of mutualistic pairs that trade it in **either**
direction (a→b ∪ b→a, de-duplicated per pair). / 把 L0–L9 所有互利菌对汇总，
统计每种交换代谢物在**任一方向**（a→b ∪ b→a，按菌对去重）被交换的比例。

---

## 6. Supplementary analyses / 补充分析

Beyond the gradient sweep, several focused analyses back individual
supplementary figures. Each has its own `results/` sub-directory and plotting
script. / 除梯度扫描外，以下针对性分析各自支撑一张补充图，均有独立的 `results/`
子目录与绘图脚本。

| Fig 图 | Analysis 分析 | Data 数据 | Takeaway 结论 |
|--------|---------------|-----------|---------------|
| **S1** | L0 substrate basis / L0 底物来源 | `results/l0_substrates/lac_L0_uptake.tsv` | All 10 *Lactobacillus* are viable at L0 on amino acids + mucin amino sugars (no fatty-acid oxidation). / 10 株在 L0 均可存活，靠氨基酸 + 黏蛋白氨基糖。 |
| **S5** | Threshold sensitivity / 阈值敏感性 | `results/{akk,lac}_vs_uhgg/` | The pectin peak / glucose crash survives every mutualism-/competition-threshold combination. / 阈值组合下结论不变。 |
| **S6** | Objective-function sensitivity / 目标函数敏感性 | `results/objective_sensitivity/` | Default **lexicographic max-min** vs **fixed-ratio** give identical mutualism fractions and ≥ 99 % per-pair agreement. / 两种共培养目标结果一致。 |
| **S7** | Controlled glucose contrast / 受控葡萄糖对照 | `results/glucose_isolation/` | Holding L1–L4 fixed, swapping pectin (L5) → glucose (L5-glc) lowers mutualism in **both** systems (Akk 10.4→6.5 %, Lac 21.7→18.4 %): a shared carbon-**quality** effect. / 同背景下果胶换葡萄糖，两系统互利率都下降，是共有的碳源**质量**效应。 |
| **S8** | GPR coverage / 基因支持率 | `results/gpr_analysis/gpr_coverage.tsv` | Cross-feeding gene support (~50–60 %) tracks transporter annotation, not a gap-filling artefact. / 交叉喂养基因支持率匹配转运体注释，非补缺假象。 |

**Re-running the objective-function & glucose analyses** (needs the fast-mic
binary) / 重跑目标函数与葡萄糖分析（需 fast-mic 二进制）：

```bash
# Objective sensitivity (Fig S6): L5/L6 under two co-culture objectives
fast-mic --group1 <probiotic_dir> --group2 test/UHGG/final_gapseq_xml \
  --medium-file <gradient_L5_pectin.csv> --summary -o results/objective_sensitivity/L5_pectin_lexicographic.tsv
fast-mic ... --fixed-ratio -o results/objective_sensitivity/L5_pectin_fixed_ratio.tsv   # repeat for L6

# Controlled glucose contrast (Fig S7): same L1–L4 background, pectin vs glucose
fast-mic --group1 test/akk/... --group2 test/UHGG/final_gapseq_xml \
  --medium-file <gradient_L5glc.csv> --summary -o results/glucose_isolation/akk_L5glc.tsv   # repeat akk/lac × L5/L5glc/L6
```

---

## 7. Data dictionary / 数据字典

Each level produces a summary `L*.tsv` (one row per species pair). With
`--full-tsv`, a richer `L*.full.tsv` adds per-metabolite cross-feeding columns.
每一级生成汇总表 `L*.tsv`（每行一个菌对）；加 `--full-tsv` 时额外生成
`L*.full.tsv`，含逐代谢物交叉喂养信息。（`lac_vs_uhgg/*.full.tsv` 以 `.gz` 压缩。）

| Column 列                     | Meaning 含义                                                       |
|-------------------------------|--------------------------------------------------------------------|
| `species_a` / `species_b`     | 菌对的两个成员 / the two members of the pair                       |
| `growth_a_alone` / `_b_alone` | 单培养生长速率 / monoculture growth rate (h⁻¹)                    |
| `growth_a_co` / `_b_co`       | 共培养生长速率 / co-culture growth rate (h⁻¹)                      |
| `benefit_a` / `benefit_b`     | 相对收益 / relative benefit (co vs mono)                           |
| `interaction_type`            | 互作类型 / mutualism, competition, commensalism, parasitism, …     |
| `competition_intensity`       | 资源重叠强度 / shared-uptake overlap                               |
| `n_exchanged_metabolites`     | 交换代谢物数 / number of cross-fed metabolites                     |
| `gene_supported_fraction`     | 有基因证据的交叉喂养通量占比 / gene-supported cross-feeding flux    |
| `a_to_b_metabolites` ⁺        | a→b 交换的代谢物/通量/供受体基因 / metabolites, fluxes, donor/receiver genes |
| `competed_metabolites` ⁺      | 双方竞争的代谢物 / metabolites competed for                        |

⁺ `*.full.tsv` only / 仅 `*.full.tsv` 含有。

A viable pair = both members grow in monoculture (> 1×10⁻⁴ h⁻¹); mutualism % is
computed over viable pairs. / 可行菌对 = 两成员单培养均可生长（> 1×10⁻⁴ h⁻¹）；
互利比例在可行菌对上统计。

---

## 8. Figure index / 图索引

| Figure 图 | Script 脚本                              | Topic 主题                                      |
|-----------|------------------------------------------|-------------------------------------------------|
| Fig 1     | `plot_fig1_benchmark.R`                  | fast-mic vs COBRApy 性能基准 / benchmark (A accuracy, B runtime, C speedup, D memory) |
| Fig 2     | `plot_fig2_phylo_trees.R`                | 益生菌系统发育树 / probiotic phylogenies        |
| Fig 3     | `plot_fig3_gradient_overview.R`          | 梯度上的互作总览 / interactions across gradient |
| Fig 4     | `plot_fig4_mechanisms.R`                 | 葡萄糖崩溃的机制 / glucose-crash mechanisms     |
| Fig 5     | `plot_fig5_crossfeed_sankey.R`           | 交叉喂养货币桑基图 / cross-feeding currencies   |
| Fig 6     | `plot_fig6_strain_heatmap.R`             | 菌株级合作异质性 / strain-level heterogeneity   |
| Fig S1    | `plot_figS1_l0_substrates.R`             | L0 底物来源 / L0 substrate basis of viability   |
| Fig S2    | `plot_figS2_interaction_composition.R`   | 六类互作组成 / six-category composition         |
| Fig S3    | `plot_figS3_benefit_landscape.R`         | 收益分布景观 / benefit landscape (β_A vs β_B)   |
| Fig S4    | `plot_figS4_crossfeed_landscape.R`       | 交叉喂养代谢物景观 / cross-feeding landscape    |
| Fig S5    | `plot_figS5_threshold_sensitivity.R`     | 阈值敏感性 / threshold robustness               |
| Fig S6    | `plot_figS6_objective_sensitivity.R`     | 目标函数敏感性 / objective-function sensitivity |
| Fig S7    | `plot_figS7_glucose_isolation.R`         | 受控葡萄糖对照 / controlled glucose contrast    |
| Fig S8    | `plot_figS8_gpr_coverage.R`              | 交叉喂养基因支持率 / cross-feeding GPR coverage |

> Supplementary figures are numbered in order of first citation in the
> manuscript. / 补充图按正文首次引用顺序编号。

---

## 9. Citation / 引用

If you use these data or scripts, please cite the fast-mic manuscript.
如使用本数据或脚本，请引用 fast-mic 论文。
- Engine 引擎: <https://github.com/Zhangjyfree/fast-mic.git>
- Reproduction 复现: <https://github.com/Zhangjyfree/fast-mic-manuscript.git>
