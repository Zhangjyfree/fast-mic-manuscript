#!/usr/bin/env Rscript
# Figure 5 — Predicted cross-feeding currencies of the two probiotics
#
# (A) Sankey/alluvial flow  probiotic -> exported metabolite -> biological function.
#     Flow width = % of mutualistic pairs (pooled over the L0-L9 gradient) that
#     exchange the metabolite, taken from the pooled cross-feeding table so the
#     numbers match the manuscript. Two genome-encoded currency axes are explicit:
#     Akkermansia exports organic acids & 1,2-propanediol; Lactobacillus exports
#     purine nucleosides & pentoses. * = independent experimental support.
# (B) For key currencies, the % of viable pairs exchanging the metabolite through a
#     mutualistic interaction at L5 (pectin) vs L6 (free glucose): the glucose crash.
#
# Usage: Rscript scripts/plot_fig5_crossfeed_sankey.R [--results results]

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggalluvial); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
RESULTS <- get_arg("--results", "results")
OUTDIR  <- get_arg("--outdir",  file.path(RESULTS, "figures_paper"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
THRESH <- 1e-4
COL_AKK <- "#762a83"; COL_LAC <- "#1b7837"
# ── currency catalogue (manuscript) ──────────────────────────────────────────
CUR <- tibble::tribble(
  ~probiotic,      ~cpd,        ~name,          ~func,                                ~ref,
  "Akkermansia",   "cpd00116",  "Methanol",     "Deoxy-sugar\nproducts",       FALSE,
  "Akkermansia",   "cpd00118",  "Putrescine",   "Deoxy-sugar\nproducts",       FALSE,
  "Akkermansia",   "cpd01861",  "(R)-1,2-PD",   "Deoxy-sugar\nproducts",       TRUE,
  "Akkermansia",   "cpd00036",  "Succinate",    "SCFA &\npropionate", TRUE,
  "Akkermansia",   "cpd00047",  "Formate",      "SCFA &\npropionate", FALSE,
  "Akkermansia",   "cpd00141",  "Propionate",   "SCFA &\npropionate", FALSE,
  "Akkermansia",   "cpd00106",  "Fumarate",     "SCFA &\npropionate", FALSE,
  "Akkermansia",   "cpd00130",  "L-Malate",     "SCFA &\npropionate", FALSE,
  "Lactobacillus", "cpd00246",  "Inosine",      "Purine\nsalvage",         TRUE,
  "Lactobacillus", "cpd01217",  "Xanthosine",   "Purine\nsalvage",         TRUE,
  "Lactobacillus", "cpd00311",  "Guanosine",    "Purine\nsalvage",         TRUE,
  "Lactobacillus", "cpd00367",  "Cytidine",     "Purine\nsalvage",         FALSE,
  "Lactobacillus", "cpd00105",  "D-Ribose",     "Pentose\nutilisation",               FALSE,
  "Lactobacillus", "cpd00130",  "L-Malate",     "SCFA &\npropionate", FALSE
)

# ── currency frequencies from the pooled table (matches manuscript, all genera) ──
cf <- suppressMessages(read_tsv(file.path(OUTDIR, "crossfeed_landscape_table.tsv"),
                                show_col_types = FALSE, progress = FALSE))
syslab <- function(pro) if (pro == "Akkermansia") "Akkermansia" else "Lactobacillus × UHGG"
CUR$freq <- mapply(function(pro, cpd) {
  r <- cf[grepl(syslab(pro), cf$sys_label, fixed = TRUE) & cf$cpd == cpd, ]
  if (nrow(r) > 0) 100 * r$frac[1] else 0
}, CUR$probiotic, CUR$cpd)

dat <- CUR |> filter(freq >= 5) |>
  mutate(met_lab = ifelse(cpd == "cpd00130", "L-Malate  (shared)",
                          sprintf("%s  %.0f%%%s", name, round(freq), ifelse(ref, " *", ""))),
         probiotic = factor(probiotic, levels = c("Akkermansia","Lactobacillus")),
         func = factor(func, levels = unique(CUR$func)))
met_lab_levels <- dat |> arrange(probiotic, desc(freq)) |> pull(met_lab) |> unique()
dat <- dat |> mutate(met_lab = factor(met_lab, levels = met_lab_levels))

FIG_W <- 232/25.4; FIG_H <- 112/25.4
pA <- ggplot(dat, aes(axis1 = probiotic, axis2 = met_lab, axis3 = func, y = freq)) +
  geom_alluvium(aes(fill = probiotic), width = 0.22, alpha = 0.48,
                curve_type = "sigmoid", colour = NA) +
  geom_stratum(aes(fill = probiotic), width = 0.22, linewidth = 0.6, colour = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.05, lineheight = 0.9, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Akkermansia" = COL_AKK, "Lactobacillus" = COL_LAC), name = NULL) +
  scale_x_discrete(limits = c("Probiotic", "Exported metabolite", "Biological function / fate"),
                   expand = expansion(mult = c(0.13, 0.16)), position = "top") +
  labs(title = "Predicted cross-feeding currencies of the two probiotics",
       subtitle = "Flow width = % of mutualistic pairs (pooled over the L0\u2013L9 gradient) exchanging the metabolite      * = independent experimental support",
       y = NULL) +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5, colour = "grey15",
                                  margin = margin(b = 1)),
        plot.subtitle = element_text(size = 5.8, colour = "grey40", hjust = 0.5, margin = margin(b = 7)),
        axis.text.x = element_text(face = "bold", size = 6.8, colour = "grey30"),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.position = "none", plot.margin = margin(8, 18, 4, 14))

# (Panel B removed: no suitable place in the draft)


ggsave(file.path(OUTDIR, "fig5_crossfeed_sankey.tiff"), pA, width = FIG_W, height = FIG_H,
       dpi = 300, compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "fig5_crossfeed_sankey.pdf"), pA, width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "fig5_crossfeed_sankey.png"), pA, width = FIG_W, height = FIG_H, dpi = 300)
cat("Saved fig5_crossfeed_sankey (currency -> function sankey) to", OUTDIR, "\n")
print(dat |> select(probiotic, name, freq, func))
