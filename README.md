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
figures (Figures 1–6 and Supplementary Figures S1–S5). To re-run the
simulations you also need the compiled `fast-mic` binary from the engine repo.

**中文** —— 本仓库收录 fast-mic 的原始模拟输出、基因组/模型输入，以及复现正文图
（Figure 1–6）和补充图（Figure S1–S5）所需的全部绘图脚本。若要重跑模拟，还需
从引擎仓库编译得到 `fast-mic` 可执行文件。

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
│   └── plot_figS1–S5_*.R               # 补充图 S1–S5 / Supplementary figures
├── results/                            # fast-mic 输出 + 成图 / outputs + figures
│   ├── akk_vs_uhgg/                    # Akkermansia × UHGG, L0–L9
│   ├── lac_vs_uhgg/                    # Lactobacillus × UHGG, L0–L9
│   ├── objective_sensitivity/          # fixed-ratio vs lexicographic (Fig S5)
│   └── figures_paper/                  # 最终 PDF + 汇总表 / final figures + tables
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
`results/` and writing PDF + PNG + TIFF into `results/figures_paper/`.

**中文** —— `results/` 中已附带模拟输出，可直接重绘所有图。每个 R 脚本默认从
`results/` 读取，并把 PDF + PNG + TIFF 输出到 `results/figures_paper/`。

```bash
# 正文图 / Main figures
Rscript scripts/plot_fig1_benchmark.R \
    benchmark/benchmark_results/thread_scaling/thread_scaling_results.tsv \
    benchmark/benchmark_results/correctness/scatter.tsv \
    benchmark/benchmark_results/correctness/stats.tsv      # ↑ benchmark data lives in the engine repo
Rscript scripts/plot_fig2_phylo_trees.R
Rscript scripts/plot_fig3_gradient_overview.R
Rscript scripts/plot_fig4_mechanisms.R
Rscript scripts/plot_fig5_crossfeed_sankey.R              # 需先生成下方汇总表 / needs the table below
Rscript scripts/plot_fig6_strain_heatmap.R

# 补充图 / Supplementary figures
Rscript scripts/plot_figS1_interaction_composition.R
Rscript scripts/plot_figS2_benefit_landscape.R
Rscript scripts/plot_figS3_crossfeed_landscape.R
Rscript scripts/plot_figS4_threshold_sensitivity.R
Rscript scripts/plot_figS5_objective_sensitivity.R
```

**Cross-feeding summary table / 交叉喂养汇总表** (used by Fig 5 & S3):

```bash
python3 scripts/make_crossfeed_table.py
# → results/figures_paper/crossfeed_landscape_table.tsv
```

It pools every mutualistic pair across L0–L9 and counts, for each exchanged
metabolite, the fraction of mutualistic pairs that trade it in **either**
direction (a→b ∪ b→a, de-duplicated per pair). / 把 L0–L9 所有互利菌对汇总，
统计每种交换代谢物在**任一方向**（a→b ∪ b→a，按菌对去重）被交换的比例。

---

## 6. Objective-function sensitivity / 目标函数敏感性 (Fig S5)

`results/objective_sensitivity/` holds *L. gasseri* × UHGG results at L5 and L6
under two co-culture objectives — the default **lexicographic max-min** and an
alternative **fixed-ratio** (`fast-mic --fixed-ratio`, which pins the co-culture
growth ratio to the monoculture ratio). The two objectives give identical
mutualism fractions and ≥ 99 % per-pair agreement, showing the conclusions are
robust to the allocation rule (Figure S5). / 此目录存放 *L. gasseri* × UHGG 在
L5/L6 下、两种共培养目标（默认字典序最大最小 vs `--fixed-ratio` 固定比例）的结果；
两者给出完全一致的互利比例与 ≥99% 的逐对一致率，证明结论不依赖分配规则（图 S5）。

```bash
# regenerate the four tables (needs the fast-mic binary) / 重新生成四张表（需 fast-mic 二进制）
fast-mic --group1 <gasseri_dir> --group2 test/UHGG/final_gapseq_xml \
  --medium-file <gradient_L5_pectin.csv> --summary -o L5_lexicographic.tsv
fast-mic ... --fixed-ratio -o L5_fixed_ratio.tsv          # repeat for L6
Rscript scripts/plot_figS5_objective_sensitivity.R
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

---

## 8. Figure index / 图索引

| Figure 图 | Script 脚本                              | Topic 主题                                      |
|-----------|------------------------------------------|-------------------------------------------------|
| Fig 1     | `plot_fig1_benchmark.R`                  | fast-mic vs COBRApy 性能基准 / benchmark        |
| Fig 2     | `plot_fig2_phylo_trees.R`                | 益生菌系统发育树 / probiotic phylogenies        |
| Fig 3     | `plot_fig3_gradient_overview.R`          | 梯度上的互作总览 / interactions across gradient |
| Fig 4     | `plot_fig4_mechanisms.R`                 | 葡萄糖崩溃的机制 / glucose-crash mechanisms     |
| Fig 5     | `plot_fig5_crossfeed_sankey.R`           | 交叉喂养货币桑基图 / cross-feeding currencies   |
| Fig 6     | `plot_fig6_strain_heatmap.R`             | 菌株级合作异质性 / strain-level heterogeneity   |
| Fig S1    | `plot_figS1_interaction_composition.R`   | 六类互作组成 / six-category composition         |
| Fig S2    | `plot_figS2_benefit_landscape.R`         | 收益分布景观 / benefit landscape (β_A vs β_B)   |
| Fig S3    | `plot_figS3_crossfeed_landscape.R`       | 交叉喂养代谢物景观 / cross-feeding landscape    |
| Fig S4    | `plot_figS4_threshold_sensitivity.R`     | 阈值敏感性 / threshold robustness               |
| Fig S5    | `plot_figS5_objective_sensitivity.R`     | 目标函数敏感性 / objective-function sensitivity |

---

## 9. Citation / 引用

If you use these data or scripts, please cite the fast-mic manuscript.
如使用本数据或脚本，请引用 fast-mic 论文。
- Engine 引擎: <https://github.com/Zhangjyfree/fast-mic.git>
- Reproduction 复现: <https://github.com/Zhangjyfree/fast-mic-manuscript.git>
