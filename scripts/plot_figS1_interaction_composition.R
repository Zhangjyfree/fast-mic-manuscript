#!/usr/bin/env Rscript
# =============================================================================
# plot_figS1_interaction_composition.R   (Supplementary Figure S1)
# -----------------------------------------------------------------------------
# Full six-category interaction-type composition across the prebiotic gradient.
# Main Figure 3 shows only mutualism and competition fractions; this panel
# resolves all six ecological categories so the reader can see that the
# remaining classes (commensalism, parasitism, amensalism, neutral) stay
# negligible across the whole gradient and every system.
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

THRESHOLD <- 1e-4
FIG_W <- 170 / 25.4
FIG_H <- 70 / 25.4
FIG_DPI <- 300

LEVELS <- tibble::tribble(
  ~level_idx, ~level_file,              ~short,
  0, "L0_base","L0", 1, "L1_inulin","L1", 2, "L2_fos","L2",
  3, "L3_gos","L3", 4, "L4_xos","L4", 5, "L5_pectin","L5",
  6, "L6_resistant_starch","L6", 7, "L7_bglucan","L7",
  8, "L8_hmo","L8", 9, "L9_mos","L9"
)

SYSTEMS <- tibble::tribble(
  ~sys_dir,        ~label,
  "akk_vs_uhgg",   "Akkermansia × Gut (UHGG)",
  "lac_vs_uhgg",   "Lactobacillus × Gut (UHGG)"
)

# Canonical order + colours for the six interaction types
ITYPES <- c("mutualism", "commensalism", "neutral",
            "parasitism", "amensalism", "competition")
ITYPE_COLOURS <- c(
  mutualism    = "#1b7837",
  commensalism = "#a6dba0",
  neutral      = "#cccccc",
  parasitism   = "#fdb863",
  amensalism   = "#e08214",
  competition  = "#762a83"
)
ITYPE_LABELS <- c(
  mutualism    = "Mutualism",
  commensalism = "Commensalism",
  neutral      = "Neutral",
  parasitism   = "Parasitism",
  amensalism   = "Amensalism",
  competition  = "Competition"
)

# ── Load: per-level composition over all six categories ──────────────────────
load_system <- function(sys_dir, label) {
  rows <- list()
  for (l in seq_len(nrow(LEVELS))) {
    p <- file.path(RESULTS, sys_dir, paste0(LEVELS$level_file[l], ".tsv"))
    if (!file.exists(p)) next
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    viable <- df |> filter(growth_a_alone > THRESHOLD, growth_b_alone > THRESHOLD)
    if (nrow(viable) == 0) next
    comp <- viable |>
      mutate(interaction_type = tolower(interaction_type)) |>
      count(interaction_type) |>
      mutate(frac = n / sum(n))
    rows[[length(rows)+1]] <- comp |>
      mutate(system = label, level_idx = LEVELS$level_idx[l])
  }
  bind_rows(rows)
}

comp <- bind_rows(map2(SYSTEMS$sys_dir, SYSTEMS$label, load_system)) |>
  mutate(interaction_type = factor(interaction_type, levels = ITYPES),
         system = factor(system, levels = SYSTEMS$label)) |>
  left_join(LEVELS, by = "level_idx") |>
  mutate(short = factor(short, levels = LEVELS$short))

# ── Theme ────────────────────────────────────────────────────────────────────
theme_gm <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title      = element_text(face = "bold", size = 7, hjust = 0),
      strip.text      = element_text(face = "bold", size = 6.5),
      strip.background = element_blank(),
      axis.title      = element_text(size = 6.5),
      axis.text       = element_text(size = 5.5, colour = "black"),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.line       = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.25, colour = "black"),
      legend.position = "right",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 6),
      legend.key.size = unit(0.3, "cm"),
      plot.margin     = margin(2, 4, 2, 2)
    )
}

glucose_vline <- geom_vline(xintercept = 5.5, linetype = "dashed",
                            colour = "#CC3311", linewidth = 0.3)

figS1 <- ggplot(comp, aes(level_idx, frac, fill = interaction_type)) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
  glucose_vline +
  facet_wrap(~ system, nrow = 1) +
  scale_fill_manual(values = ITYPE_COLOURS, labels = ITYPE_LABELS,
                    breaks = ITYPES) +
  scale_x_continuous(breaks = LEVELS$level_idx, labels = LEVELS$short,
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Prebiotic gradient (L0–L9)",
       y = "Fraction of viable pairs",
       title = "Full interaction-type composition across the prebiotic gradient") +
  theme_gm()

ggsave(file.path(OUTDIR, "figS1_interaction_composition.tiff"), figS1,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "figS1_interaction_composition.pdf"), figS1,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "figS1_interaction_composition.png"), figS1,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved figS1_interaction_composition (TIFF/PDF/PNG) to", OUTDIR, "\n")
