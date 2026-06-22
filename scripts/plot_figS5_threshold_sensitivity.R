#!/usr/bin/env Rscript
# =============================================================================
# plot_figS5_threshold_sensitivity.R   (Supplementary Figure S5)
# -----------------------------------------------------------------------------
# Robustness of the central result (inverted-U mutualism peaking at L5,
# collapsing at L6) to the two classification cut-offs:
#   A  viability threshold  mu_min  in {1e-5, 1e-4, 1e-3}
#   B  benefit threshold    eps_beta in {1e-4, 1e-3, 1e-2}
# Mutualism is re-derived directly from benefit_a / benefit_b, so the figure
# shows the qualitative pattern is preserved across an order of magnitude in
# either parameter.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
RESULTS <- get_arg("--results", "results")
OUTDIR  <- get_arg("--outdir",  file.path(RESULTS, "figures_paper"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FIG_W <- 170 / 25.4
FIG_H <- 80 / 25.4
FIG_DPI <- 300

LEVELS <- tibble::tribble(
  ~level_idx, ~level_file,              ~short,
  0, "L0_base","L0", 1, "L1_inulin","L1", 2, "L2_fos","L2",
  3, "L3_gos","L3", 4, "L4_xos","L4", 5, "L5_pectin","L5",
  6, "L6_resistant_starch","L6", 7, "L7_bglucan","L7",
  8, "L8_hmo","L8", 9, "L9_mos","L9"
)

SYSTEMS <- tibble::tribble(
  ~sys_dir,        ~label,                            ~colour,
  "akk_vs_uhgg",   "Akkermansia × Gut (UHGG)",        "#762a83",
  "lac_vs_uhgg",   "Lactobacillus × Gut (UHGG)",      "#1b7837"
)
sys_colours <- setNames(SYSTEMS$colour, SYSTEMS$label)

VIAB_SET <- c(1e-5, 1e-4, 1e-3)
EPS_SET  <- c(1e-4, 1e-3, 1e-2)
DEFAULT_VIAB <- 1e-4
DEFAULT_EPS  <- 1e-3

# ── Load raw benefits once per system × level ────────────────────────────────
load_raw <- function(sys_dir, label) {
  rows <- list()
  for (l in seq_len(nrow(LEVELS))) {
    p <- file.path(RESULTS, sys_dir, paste0(LEVELS$level_file[l], ".tsv"))
    if (!file.exists(p)) next
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    rows[[length(rows)+1]] <- df |>
      transmute(system = label, level_idx = LEVELS$level_idx[l],
                ga = growth_a_alone, gb = growth_b_alone,
                ba = benefit_a, bb = benefit_b)
  }
  bind_rows(rows)
}

raw <- bind_rows(map2(SYSTEMS$sys_dir, SYSTEMS$label, load_raw)) |>
  mutate(system = factor(system, levels = SYSTEMS$label))

mut_frac <- function(d, viab, eps) {
  d |>
    filter(ga > viab, gb > viab) |>
    group_by(system, level_idx) |>
    summarise(frac = mean(ba > eps & bb > eps), .groups = "drop")
}

# A: vary viability threshold (eps fixed)
dfA <- map_dfr(VIAB_SET, function(v)
  mut_frac(raw, v, DEFAULT_EPS) |> mutate(param = v)) |>
  mutate(param = factor(sprintf("%g", param), levels = sprintf("%g", VIAB_SET)))

# B: vary benefit threshold (viab fixed)
dfB <- map_dfr(EPS_SET, function(e)
  mut_frac(raw, DEFAULT_VIAB, e) |> mutate(param = e)) |>
  mutate(param = factor(sprintf("%g", param), levels = sprintf("%g", EPS_SET)))

theme_gm <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title       = element_text(face = "bold", size = 7, hjust = 0),
      axis.title       = element_text(size = 6.5),
      axis.text        = element_text(size = 5.5, colour = "black"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.15),
      legend.position  = c(0.5, 0.93),
      legend.direction = "horizontal",
      legend.title     = element_text(size = 5.5),
      legend.text      = element_text(size = 5),
      legend.key.size  = unit(0.22, "cm"),
      legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
      plot.margin      = margin(2, 4, 2, 2)
    )
}

x_sc <- scale_x_continuous(breaks = LEVELS$level_idx, labels = LEVELS$short,
                           expand = expansion(add = 0.3))
glucose_vline <- geom_vline(xintercept = 5.5, linetype = "dashed",
                            colour = "#CC3311", linewidth = 0.3)

mk_panel <- function(df, lty_name, title) {
  ggplot(df, aes(level_idx, frac, colour = system, linetype = param)) +
    glucose_vline +
    geom_line(linewidth = 0.45) + geom_point(size = 0.7) +
    scale_colour_manual(values = sys_colours, guide = "none") +
    scale_linetype_manual(values = c("dotted", "solid", "dashed"),
                          name = lty_name) +
    x_sc +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    labs(title = title, x = "Prebiotic gradient (L0–L9)",
         y = "Mutualism (% of viable pairs)") +
    theme_gm()
}

pA <- mk_panel(dfA, expression(mu[min]),
               "A  Viability-threshold sensitivity")
pB <- mk_panel(dfB, expression(epsilon[beta]),
               "B  Benefit-threshold sensitivity")

figS5 <- pA + pB + plot_layout(nrow = 1) &
  theme(plot.tag = element_text(face = "bold", size = 9, family = "Arial"))

ggsave(file.path(OUTDIR, "figS5_threshold_sensitivity.tiff"), figS5,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "figS5_threshold_sensitivity.pdf"), figS5,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "figS5_threshold_sensitivity.png"), figS5,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved figS5_threshold_sensitivity (TIFF/PDF/PNG) to", OUTDIR, "\n")
