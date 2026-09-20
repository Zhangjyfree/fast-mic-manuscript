#!/usr/bin/env Rscript
# Fig 5 (combined): predicted cross-feeding currencies of the two probiotics.
#   A  Sankey/alluvial  probiotic -> exported metabolite -> biological fate (original)
#   B  Condition-resolved prevalence: top exported metabolites x L0-L9 (Akk & Lac)
#   C  Metabolite-class composition of each probiotic's cross-feeding currencies
# Ambiguous universal by-products (ethanol, pyruvate, L/D-lactate) are set aside
# (as in the manuscript) for B and C. Sankey frequencies are the pooled % of
# mutualistic pairs, computed from the same extracted table.
#
# Input: results/fig5/fig5_prevalence_{lac,akk}.tsv  (extract_fig5_crossfeed.py)
# Usage: Rscript scripts/plot_fig5_v2.R [output_prefix]

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggalluvial); library(patchwork); library(scales)
})

## shared figure conventions (BMC/Microbiome): Arial, 170 mm geometry,
## italic species names.  See scripts/figure_theme.R
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  cand <- c(if (length(f)) file.path(dirname(f[1]), "figure_theme.R"),
            "scripts/figure_theme.R", "figure_theme.R")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("figure_theme.R not found")
  source(hit[1], local = FALSE)
})


args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[1] else "results/figures_paper/fig5_v2_crossfeed"
FM <- "."   # repo root

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 200 / 25.4      # within the 225 mm limit
LEVELS <- c("L0","L1","L2","L3","L4","L5","L6","L7","L8","L9")
AMBIG  <- c("cpd00363","cpd00020","cpd00159","cpd00221")
TOPN   <- 12
COL_AKK <- "#762A83"; COL_LAC <- "#1B7837"
CLASS_COL <- c("organic acid"="#D6604D","nucleoside"="#4393C3","sugar"="#5AAE61",
               "amino acid"="#F4A582","alcohol/diol"="#9970AB","other"="#BBBBBB")
SYS_LAB <- c(akk="Akkermansia", lac="Lactobacillus")
GENERA  <- unname(SYS_LAB)   # stratum labels that are genus names

base_theme <- theme_bw(base_size = 9, base_family = BMC_FONT) +
  theme(panel.grid = element_blank(), legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 5.8), legend.title = element_text(size = 6.3),
        axis.text = md(size = 6), axis.text.x = md(size = 6),
        axis.text.y = md(size = 6), axis.title = element_text(size = 7),
        plot.title = md(size = 8, face = "bold"), strip.text = md(size = 7))

# ModelSEED ships a few compounds under a terse abbreviation rather than a name
# (MTTL, PAN, XAN ...). Spell them out so the panels match Table S11, where the
# reviewer asked for these abbreviations to be defined.
NAME_FIX <- c("cpd00324" = "Methanethiol",   "cpd00644" = "Pantothenate",
              "cpd00309" = "Xanthine",       "cpd00208" = "Lactose",
              "cpd00222" = "D-Gluconate",    "cpd00276" = "D-Glucosamine",
              "cpd00359" = "Indole")

read_sys <- function(s) read_tsv(sprintf("%s/results/fig5/fig5_prevalence_%s.tsv", FM, s),
                                 show_col_types = FALSE) |>
  mutate(system = SYS_LAB[s], sys = s, prevalence = as.numeric(prevalence),
         name = ifelse(cpd %in% names(NAME_FIX), NAME_FIX[cpd], name),
         level = factor(level, levels = LEVELS))
raw <- bind_rows(read_sys("akk"), read_sys("lac"))
d   <- raw |> filter(!cpd %in% AMBIG)

# ── pooled frequency (% of mutualistic pairs over the whole gradient) per cpd ──
pooled <- raw |> group_by(sys, cpd, name) |>
  summarise(exp = sum(n_export), .groups = "drop") |>
  left_join(raw |> distinct(sys, level, n_mut_pairs) |>
              group_by(sys) |> summarise(totmut = sum(n_mut_pairs), .groups = "drop"),
            by = "sys") |>
  mutate(freq = 100 * exp / totmut)

# ── Panel A: Sankey (currency catalogue from manuscript; freq from pooled) ────
CUR <- tibble::tribble(
  ~probiotic,      ~sys,  ~cpd,        ~name,         ~func,
  "Akkermansia",   "akk", "cpd00116",  "Methanol",    "Deoxy-sugar\nproducts",
  "Akkermansia",   "akk", "cpd00118",  "Putrescine",  "Deoxy-sugar\nproducts",
  "Akkermansia",   "akk", "cpd01861",  "(R)-1,2-PD*", "Deoxy-sugar\nproducts",
  "Akkermansia",   "akk", "cpd00036",  "Succinate*",  "SCFA &\npropionate",
  "Akkermansia",   "akk", "cpd00047",  "Formate",     "SCFA &\npropionate",
  "Akkermansia",   "akk", "cpd00141",  "Propionate",  "SCFA &\npropionate",
  "Akkermansia",   "akk", "cpd00130",  "L-Malate",    "SCFA &\npropionate",
  "Lactobacillus", "lac", "cpd00246",  "Inosine*",    "Purine\nsalvage",
  "Lactobacillus", "lac", "cpd01217",  "Xanthosine*", "Purine\nsalvage",
  "Lactobacillus", "lac", "cpd00311",  "Guanosine*",  "Purine\nsalvage",
  "Lactobacillus", "lac", "cpd00367",  "Cytidine",    "Purine\nsalvage",
  "Lactobacillus", "lac", "cpd00105",  "D-Ribose",    "Pentose\nutilisation"
) |>
  left_join(pooled |> select(sys, cpd, freq), by = c("sys","cpd")) |>
  mutate(freq = ifelse(is.na(freq), 0, freq)) |>
  filter(freq >= 5) |>
  mutate(met_lab = sprintf("%s  %.0f%%", name, round(freq)),
         probiotic = factor(probiotic, levels = c("Akkermansia","Lactobacillus")),
         func = factor(func, levels = unique(func)))
CUR <- CUR |> mutate(met_lab = factor(met_lab, levels = CUR |> arrange(probiotic, desc(freq)) |> pull(met_lab)))

pA <- ggplot(CUR, aes(axis1 = probiotic, axis2 = met_lab, axis3 = func, y = freq)) +
  geom_alluvium(aes(fill = probiotic), width = 0.36, alpha = 0.5, curve_type = "sigmoid", colour = NA) +
  geom_stratum(aes(fill = probiotic), width = 0.36, linewidth = 0.5, colour = "white") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum),
                fontface = after_stat(ifelse(as.character(stratum) %in% GENERA,
                                             "bold.italic", "bold"))),
            size = 1.45, lineheight = 0.85, colour = "white") +
  scale_discrete_identity(aesthetics = "fontface") +
  scale_fill_manual(values = c("Akkermansia"=COL_AKK,"Lactobacillus"=COL_LAC), guide = "none") +
  scale_x_discrete(limits = c("Probiotic","Exported metabolite","Biological fate"),
                   expand = expansion(mult = c(0.10, 0.12)), position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.04))) +   # room for the last stratum label
  labs(title = "Predicted cross-feeding currencies (pooled L0–L9; * = experimental support)", y = NULL) +
  theme_void(base_size = 7, base_family = BMC_FONT) +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5, margin = margin(b = 3)),
        axis.text.x = element_text(face = "bold", size = 6, colour = "grey30"),
        plot.margin = margin(4, 8, 8, 8))

# ── Panel B: condition-resolved prevalence heatmaps ──────────────────────────
heat_panel <- function(sys_name, fillcol) {
  dd <- d |> filter(system == sys_name)
  top <- dd |> group_by(cpd, name) |> summarise(m = max(prevalence), .groups = "drop") |>
    slice_max(m, n = TOPN) |> arrange(m)
  dd <- dd |> filter(cpd %in% top$cpd) |> mutate(name = factor(name, levels = top$name))
  ggplot(dd, aes(level, name, fill = prevalence)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient(low = "grey93", high = fillcol, name = "Prevalence\n(%)",
                        guide = guide_colourbar(barwidth = unit(0.25,"cm"), barheight = unit(1.8,"cm"))) +
    labs(x = NULL, y = NULL, title = sp_md(sys_name)) +
    base_theme + theme(legend.position = "right")
}
pB_akk <- heat_panel("Akkermansia", COL_AKK)
pB_lac <- heat_panel("Lactobacillus", COL_LAC)

# ── Panel C: class composition ───────────────────────────────────────────────
comp <- d |> group_by(system, class) |> summarise(w = sum(prevalence), .groups = "drop") |>
  group_by(system) |> mutate(frac = w / sum(w)) |> ungroup() |>
  mutate(class = factor(class, levels = names(CLASS_COL)))
pC <- ggplot(comp, aes(system, frac, fill = class)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(frac >= 0.06, percent(frac, accuracy = 1), "")),
            position = position_stack(vjust = 0.5), size = 1.9, colour = "white") +
  scale_fill_manual(values = CLASS_COL, name = "Metabolite class") +
  scale_x_discrete(labels = sp_md) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Share of currencies", title = "Currency class composition") +
  base_theme + theme(legend.position = "right")

# ── assemble: A on top, then (B_akk / B_lac) | C ─────────────────────────────
combined <- pA / ((pB_akk / pB_lac) | pC) +
  plot_layout(heights = c(1.45, 1.7)) +   # panel A needs height for the thin strata
  plot_annotation(tag_levels = list(c("A","B","","C")))

for (ext in c("pdf","png","tiff")) {
  path <- paste0(out_prefix, ".", ext)
  if (ext == "tiff") ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in",
                            dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "png") ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in", dpi = 300)
  else ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in", device = cairo_pdf)
  cat("Saved:", path, "\n")
}
