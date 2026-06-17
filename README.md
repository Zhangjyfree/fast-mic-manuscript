# fast-mic-manuscript

**Reproducibility repository for the fast-mic manuscript** — pairwise metabolic
interactions between probiotics (*Akkermansia*, *Lactobacillus*-group) and the
human gut microbiome (UHGG) across a 10-step prebiotic gradient.

**fast-mic 论文的复现仓库** —— 益生菌（*Akkermansia*、*Lactobacillus* 类群）
与人体肠道菌群（UHGG）之间的成对代谢互作，跨越 10 级益生元梯度（L0–L9）。

---

## 1. Overview / 概览

**EN** — This repo bundles the raw fast-mic simulation outputs, the genome /
model inputs, and every plotting script needed to regenerate the manuscript
figures (Figures 1–6 and Supplementary Figures S1–S4). The fast-mic engine
itself lives in the parent repository; here we ship only the **data + scripts**
required for reproduction.

**中文** —— 本仓库收录了 fast-mic 的原始模拟输出、基因组/模型输入，以及复现
正文图（Figure 1–6）和补充图（Figure S1–S4）所需的全部绘图脚本。fast-mic
计算引擎本体位于上级仓库；此处仅提供复现所需的**数据 + 脚本**。

---

## 2. Directory layout / 目录结构

```
fast-mic-manuscript/
├── README.md
├── scripts/                          # 全部脚本 / all scripts
│   ├── run_gradient.sh               # 跑 fast-mic 梯度模拟 / run the gradient simulation
│   ├── make_crossfeed_table.py       # 汇总交叉喂养表 / pool the cross-feeding table
│   ├── plot_fig1_benchmark.R         # 图1 性能基准 / Fig 1 benchmark
│   ├── plot_fig2_phylo_trees.R       # 图2 系统发育树 / Fig 2 phylogeny
│   ├── plot_fig3_gradient_overview.R # 图3 梯度总览 / Fig 3 gradient overview
│   ├── plot_fig4_mechanisms.R        # 图4 机制 / Fig 4 mechanisms
│   ├── plot_fig5_crossfeed_sankey.R  # 图5 交叉喂养桑基图 / Fig 5 Sankey
│   ├── plot_fig6_strain_heatmap.R    # 图6 菌株热图 / Fig 6 strain heatmap
│   └── plot_figS1–S4_*.R             # 补充图 / Supplementary figures
├── results/                          # fast-mic 输出 + 成图 / outputs + figures
│   ├── akk_vs_uhgg/                  # Akkermansia × UHGG, L0–L9
│   ├── lac_vs_uhgg/                  # Lactobacillus × UHGG, L0–L9
│   └── figures_paper/                # 最终 PDF/PNG/TIFF + 汇总表 / final figures + tables
└── test/                             # 输入基因组与模型 / input genomes & models
    ├── akk/                          # Akkermansia 基因组 / gapseq 模型 / 树
    ├── lac/                          # Lactobacillus 基因组 / gapseq 模型 / 树
    └── UHGG/                         # 肠道群落模型 + 元数据 / community models + metadata
```

**Two systems / 两个系统**

| Key 键        | Probiotic 益生菌                     | Community 群落 |
|---------------|--------------------------------------|----------------|
| `akk_vs_uhgg` | *Akkermansia* (6 strains, 3 species) | UHGG (gut)     |
| `lac_vs_uhgg` | *Lactobacillus*-group (10 strains)   | UHGG (gut)     |

---

## 3. Prebiotic gradient / 益生元梯度 (L0–L9)

**EN** — Ten strictly nested, cumulative media. Each level adds the *in vivo*
hydrolysis products of one more prebiotic on top of all previous levels. **L6 is
the first level to release free glucose** ("glucose crash").

**中文** —— 十级严格嵌套、累加式培养基。每一级在前面所有级别的基础上，再加入
一种益生元的体内水解产物。**L6 是首个释放游离葡萄糖的级别**（“葡萄糖崩溃”）。

| Level | File 文件             | Added carbon 新增碳源 / 益生元                              |
|-------|-----------------------|-----------------------------------------------------------|
| L0    | `L0_base`             | 仅基础营养 + 黏蛋白聚糖 / housekeeping + mucin glycans     |
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

**fast-mic binary** — build from the parent repo (`cargo build --release`);
`run_gradient.sh` expects it at `target/release/fast-mic`.
从上级仓库编译（`cargo build --release`），脚本默认在 `target/release/fast-mic` 调用。

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
third-party packages. 仅用标准库，无需额外安装。

---

## 5. Reproduce the figures / 复现图表

**EN** — The simulation outputs in `results/` are already provided, so you can
regenerate every figure directly. Each R script defaults to reading from
`results/` and writing PDF + PNG + TIFF into `results/figures_paper/`.

**中文** —— `results/` 中已附带模拟输出，可直接重绘所有图。每个 R 脚本默认从
`results/` 读取，并把 PDF + PNG + TIFF 输出到 `results/figures_paper/`。

```bash
# 正文图 / Main figures
Rscript scripts/plot_fig1_benchmark.R \
    benchmark/benchmark_results/thread_scaling/thread_scaling_results.tsv \
    benchmark/benchmark_results/correctness/scatter.tsv \
    benchmark/benchmark_results/correctness/stats.tsv
Rscript scripts/plot_fig2_phylo_trees.R
Rscript scripts/plot_fig3_gradient_overview.R
Rscript scripts/plot_fig4_mechanisms.R
Rscript scripts/plot_fig5_crossfeed_sankey.R          # 需先生成下方汇总表 / needs the table below
Rscript scripts/plot_fig6_strain_heatmap.R

# 补充图 / Supplementary figures
Rscript scripts/plot_figS1_threshold_sensitivity.R
Rscript scripts/plot_figS2_interaction_composition.R
Rscript scripts/plot_figS3_benefit_landscape.R
Rscript scripts/plot_figS4_crossfeed_landscape.R
```

**Cross-feeding summary table / 交叉喂养汇总表** (used by Fig 5 & S4):

```bash
python3 scripts/make_crossfeed_table.py
# → results/figures_paper/crossfeed_landscape_table.tsv
```

It pools every mutualistic pair across L0–L9 and counts, for each exchanged
metabolite, the fraction of mutualistic pairs that trade it in **either**
direction (a→b ∪ b→a, de-duplicated per pair).
它把 L0–L9 中所有互利共生的菌对汇总，统计每种交换代谢物在**任一方向**
（a→b ∪ b→a，按菌对去重）被交换的比例。

---

## 6. Re-run the simulation (optional) / 重跑模拟（可选）

**EN** — Only needed if you want to regenerate `results/*/L*.tsv` from the
genome-scale models in `test/`. Requires the compiled fast-mic binary.

**中文** —— 仅在你想从 `test/` 中的基因组尺度模型重新生成 `results/*/L*.tsv`
时才需要，依赖已编译的 fast-mic 可执行文件。

```bash
# Akkermansia × UHGG
bash scripts/run_gradient.sh \
     --group1 test/akk/akk_genomes_faa_gapseq_wdm_xml \
     --group2 test/UHGG/final_gapseq_xml \
     --threads 12 --full-tsv --out results/akk_vs_uhgg

# Lactobacillus × UHGG
bash scripts/run_gradient.sh \
     --group1 test/lac/lac_genomes_faa_gapseq_wdm_xml \
     --group2 test/UHGG/final_gapseq_xml \
     --threads 12 --full-tsv --out results/lac_vs_uhgg
```

---

## 7. Data dictionary / 数据字典

Each level produces a summary `L*.tsv` (one row per species pair). With
`--full-tsv`, a richer `L*.full.tsv` adds per-metabolite cross-feeding columns.
每一级生成汇总表 `L*.tsv`（每行一个菌对）；加 `--full-tsv` 时额外生成
`L*.full.tsv`，包含逐代谢物的交叉喂养信息。
（`lac_vs_uhgg/*.full.tsv` 以 `.gz` 压缩存储 / stored gzip-compressed.）

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

| Figure 图 | Script 脚本                            | Topic 主题                                      |
|-----------|----------------------------------------|-------------------------------------------------|
| Fig 1     | `plot_fig1_benchmark.R`                | fast-mic vs COBRApy 性能基准 / benchmark        |
| Fig 2     | `plot_fig2_phylo_trees.R`              | 益生菌系统发育树 / probiotic phylogenies        |
| Fig 3     | `plot_fig3_gradient_overview.R`        | 梯度上的互作总览 / interactions across gradient |
| Fig 4     | `plot_fig4_mechanisms.R`               | 葡萄糖崩溃的机制 / glucose-crash mechanisms     |
| Fig 5     | `plot_fig5_crossfeed_sankey.R`         | 交叉喂养货币桑基图 / cross-feeding currencies   |
| Fig 6     | `plot_fig6_strain_heatmap.R`           | 菌株级合作异质性 / strain-level heterogeneity   |
| Fig S1    | `plot_figS1_threshold_sensitivity.R`   | 阈值敏感性 / threshold sensitivity              |
| Fig S2    | `plot_figS2_interaction_composition.R` | 互作类型组成 / interaction-type composition     |
| Fig S3    | `plot_figS3_benefit_landscape.R`       | 收益分布景观 / benefit landscape                |
| Fig S4    | `plot_figS4_crossfeed_landscape.R`     | 交叉喂养代谢物景观 / cross-feeding landscape    |

---

## 9. Citation / 引用

If you use these data or scripts, please cite the fast-mic manuscript.
如使用本数据或脚本，请引用 fast-mic 论文。
fast-mic: <https://github.com/Zhangjyfree/fast-mic.git>
