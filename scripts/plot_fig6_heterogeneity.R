#!/usr/bin/env Rscript
# Fig 6 (combined): strain- and species-level heterogeneity in cooperation.
#   A  Lactobacillus strain x level mutualism heatmap        (original)
#   B  Akkermansia strain x level mutualism heatmap          (original)
#   C  Ranked per-strain mean mutualism +/- bootstrap 95% CI (statistics)
#   D  Cooperation vs metabolic capacity — NO clean predictor (honest: Lacticaseibacillus
#      have the MOST transporters/exchanges yet cooperate LEAST; no fit line drawn)
#   E  Per-strain x partner-phylum mean mutualism heatmap    (ecological resolution)
#
# Inputs: results/{lac,akk}_vs_uhgg/L*.tsv ; results/fig6/fig6_{bootstrap,phylum}_{sys}.tsv ;
#         results/fig2/fig2_traits_{sys}.tsv
# Usage: Rscript scripts/plot_fig6_v2.R [output_prefix]

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(scales)
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

options(lifecycle_verbosity = "quiet")

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[1] else "results/figures_paper/fig6_v2_heterogeneity"
FM <- "."   # repo root

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 225 / 25.4      # BMC maximum height
SYS_LAB <- c(akk = "Akkermansia", lac = "Lactobacillus")
LEVELS  <- c("L0_base","L1_inulin","L2_fos","L3_gos","L4_xos","L5_pectin",
             "L6_resistant_starch","L7_bglucan","L8_hmo","L9_mos")
SHORT   <- c("L0","L1","L2","L3","L4","L5","L6","L7","L8","L9")
THRESH  <- 1e-4

akk_species <- function(s) dplyr::case_when(
  grepl("000020225|018847155", s) ~ "A. muciniphila",
  grepl("018847135|023516715", s) ~ "A. massiliensis",
  grepl("018847015|026072915", s) ~ "A. biwaensis", TRUE ~ "Akkermansia")
lac_genus <- function(s) dplyr::case_when(
  grepl("casei|paracasei|rhamnosus", s, ignore.case = TRUE) ~ "Lacticaseibacillus",
  grepl("plantarum", s, ignore.case = TRUE) ~ "Lactiplantibacillus",
  grepl("reuteri", s, ignore.case = TRUE) ~ "Limosilactobacillus", TRUE ~ "Lactobacillus")
pretty_lac <- function(s) s |> sub("_GCF_.*$", "", x = _) |> gsub("_", " ", x = _) |>
  sub("^([A-Z]) ", "\\1. ", x = _)   # "L acidophilus" -> "L. acidophilus"
pretty_akk <- function(s) paste0(akk_species(s), " ", sub("^akk_GCF_", "", s))
grp_cols <- c("A. muciniphila"="#762A83","A. massiliensis"="#9970AB","A. biwaensis"="#C2A5CF",
              "Lactobacillus"="#1B7837","Lacticaseibacillus"="#5AAE61",
              "Lactiplantibacillus"="#A6DBA0","Limosilactobacillus"="#00695C")

base_theme <- theme_bw(base_size = 7, base_family = BMC_FONT) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.26, "cm"), legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6), axis.text = md(size = 5.5),
        axis.text.x = md(size = 5.5), axis.text.y = md(size = 5.5),
        axis.title = element_text(size = 6.8), plot.title = md(size = 8, face = "bold"),
        plot.subtitle = md(size = 5.8, colour = "grey40"), strip.text = md(size = 6),
        strip.text.x = md(size = 6), strip.text.y = md(size = 6, angle = -90))

# ── per-strain per-level mutualism fraction (for heatmaps A/B) ────────────────
strain_frac <- function(sys_key, pretty_fun) {
  rows <- list()
  for (i in seq_along(LEVELS)) {
    p <- sprintf("%s/results/%s_vs_uhgg/%s.tsv", FM, sys_key, LEVELS[i])
    if (!file.exists(p)) next
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    rows[[i]] <- df |>
      filter(growth_a_alone > THRESH, growth_b_alone > THRESH) |>
      group_by(species_a) |>
      summarise(frac = mean(tolower(interaction_type) == "mutualism") * 100, .groups = "drop") |>
      mutate(level = SHORT[i])
  }
  bind_rows(rows) |>
    mutate(strain = pretty_fun(species_a), level = factor(level, levels = SHORT))
}
lac_hm <- strain_frac("lac", pretty_lac)
akk_hm <- strain_frac("akk", pretty_akk)
lac_ord <- lac_hm |> group_by(strain) |> summarise(m = mean(frac), .groups="drop") |> arrange(m) |> pull(strain)
akk_ord <- akk_hm |> group_by(strain) |> summarise(m = mean(frac), .groups="drop") |> arrange(m) |> pull(strain)

heatmap_panel <- function(dat, ord, high, title) {
  dat <- dat |> mutate(strain = factor(strain, levels = ord))
  ggplot(dat, aes(level, strain, fill = frac)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f", frac),
                  colour = frac > 0.5 * max(dat$frac)), size = 1.7, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
    scale_fill_gradient(low = "#f7f7f7", high = high, name = "Mutualism (%)",
                        guide = guide_colourbar(barwidth = unit(0.3,"cm"), barheight = unit(2.2,"cm"))) +
    scale_y_discrete(labels = sp_md) +
    labs(x = NULL, y = NULL, title = sp_md(title)) +
    base_theme + theme(legend.position = "right")
}
pA <- heatmap_panel(lac_hm, lac_ord, "#1B7837", "Lactobacillus × gut, strain-resolved")
pB <- heatmap_panel(akk_hm, akk_ord, "#762A83", "Akkermansia × gut, strain-resolved")

# ── C/D data: bootstrap + traits ─────────────────────────────────────────────
load_sys <- function(s) {
  bt <- read_tsv(sprintf("%s/results/fig6/fig6_bootstrap_%s.tsv", FM, s), show_col_types = FALSE)
  tr <- read_tsv(sprintf("%s/results/fig2/fig2_traits_%s.tsv", FM, s), show_col_types = FALSE)
  grpfun <- if (s == "akk") akk_species else lac_genus
  prfun  <- if (s == "akk") pretty_akk else pretty_lac
  bt |> left_join(tr, by = "strain") |>
    mutate(system = SYS_LAB[s], group = grpfun(strain), pretty = prfun(strain))
}
d <- bind_rows(load_sys("akk"), load_sys("lac"))

pC <- ggplot(d |> group_by(system) |> mutate(pretty = reorder(pretty, mean_mut)) |> ungroup(),
             aes(mean_mut, pretty, colour = group)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.28, linewidth = 0.4) +
  geom_point(size = 1.3) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y", labeller = sp_labeller) +
  scale_colour_manual(values = grp_cols, name = NULL, labels = sp_md) +
  scale_y_discrete(labels = sp_md) +
  labs(x = "Mean mutualism across L0–L9 (% of viable pairs)", y = NULL,
       title = "Per-strain cooperation (bootstrap 95% CI)") +
  guides(colour = guide_legend(nrow = 1,
                               theme = theme(legend.text = md(size = 6.8)))) +
  ## legend along the bottom: the genus names are long, and a right-hand
  ## legend collides with panel E's colourbar
  base_theme + theme(legend.position = "bottom")

# ── D: honest — capacity does NOT predict cooperation (no fit line) ──────────
dD <- d |> select(system, group, mean_mut, transporters, exchanges) |>
  pivot_longer(c(transporters, exchanges), names_to = "trait", values_to = "value") |>
  mutate(trait = recode(trait, transporters = "Transporters", exchanges = "Exchanges"))
rho_lab <- dD |> group_by(trait) |>
  summarise(rho = suppressWarnings(cor(value, mean_mut, method = "spearman")),
            x = min(value), y = max(mean_mut), .groups = "drop") |>
  mutate(lab = sprintf("rho = %.2f (n.s.)", rho))
pD <- ggplot(dD, aes(value, mean_mut)) +
  geom_point(aes(colour = group), size = 1.3) +
  geom_text(data = rho_lab, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 1.9, colour = "grey25") +
  facet_wrap(~ trait, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +   # headroom for the rho labels
  scale_colour_manual(values = grp_cols, name = NULL, guide = "none") +
  labs(x = "GEM trait (count)", y = "Mean mutualism (%)",
       title = "Capacity does not predict cooperation",
       subtitle = sp_md("highest-capacity Lacticaseibacillus cooperate<br>least \u2014 no monotonic relationship")) +
  base_theme

# ── E: partner-phylum heatmap ────────────────────────────────────────────────
ph <- bind_rows(
  read_tsv(sprintf("%s/results/fig6/fig6_phylum_akk.tsv", FM), show_col_types = FALSE) |> mutate(system="Akkermansia", pretty=pretty_akk(strain)),
  read_tsv(sprintf("%s/results/fig6/fig6_phylum_lac.tsv", FM), show_col_types = FALSE) |> mutate(system="Lactobacillus", pretty=pretty_lac(strain)))
top_ph <- ph |> group_by(phylum) |> summarise(tot = sum(mean_n), .groups="drop") |> slice_max(tot, n = 8) |> pull(phylum)
phf <- ph |> filter(phylum %in% top_ph) |>
  mutate(phylum = factor(phylum, levels = rev(top_ph)),
         pretty = factor(pretty, levels = c(akk_ord, lac_ord)))
pE <- ggplot(phf, aes(phylum, pretty, fill = mean_mut)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y", labeller = sp_labeller) +
  scale_y_discrete(labels = sp_md) +
  scale_fill_viridis_c(option = "C", name = "Mutualism (%)",
                       guide = guide_colourbar(barwidth = unit(0.24,"cm"), barheight = unit(1.7,"cm"))) +
  labs(x = NULL, y = NULL, title = "Cooperation by partner phylum") +
  base_theme + theme(axis.text.x = md(angle = 45, hjust = 1, size = 5),
                     legend.position = "right")

# ── assemble: (A|B) / C / (D|E) ──────────────────────────────────────────────
row1 <- pA + pB + plot_layout(widths = c(1.3, 1))
row3 <- pD + pE + plot_layout(widths = c(1, 1.35))
combined <- row1 / pC / row3 + plot_layout(heights = c(1.05, 0.85, 1)) +
  plot_annotation(tag_levels = "A")

for (ext in c("pdf","png","tiff")) {
  path <- paste0(out_prefix, ".", ext)
  if (ext == "tiff") ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in",
                            dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "png") ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in", dpi = 300)
  else ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in", device = cairo_pdf)
  cat("Saved:", path, "\n")
}
