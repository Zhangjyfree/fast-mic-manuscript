#!/usr/bin/env Rscript
# =============================================================================
# plot_figS3_benefit_landscape.R   (Supplementary Figure S3)
# -----------------------------------------------------------------------------
# Joint distribution of the two species' relative benefits (beta_A vs beta_B)
# at representative gradient levels. The quadrants defined by epsilon_beta map
# directly onto the interaction-type classification, so this panel shows *why*
# competition dominates (both-negative quadrant is densest) and how the
# both-positive (mutualism) quadrant fills in at specialty sugars and empties
# again once glucose appears at L6.
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
EPS_BETA  <- 1e-3
FIG_W <- 170 / 25.4
FIG_H <- 150 / 25.4
FIG_DPI <- 300

# Representative levels: baseline, mutualism peak, glucose crash, saturated end
SEL <- tibble::tribble(
  ~level_file,             ~short,
  "L0_base",               "L0 Base",
  "L5_pectin",             "L5 Pectin (peak)",
  "L6_resistant_starch",   "L6 RS (+glucose)",
  "L9_mos",                "L9 MOS"
)

SYSTEMS <- tibble::tribble(
  ~sys_dir,        ~label,
  "akk_vs_uhgg",   "Akkermansia × Gut (UHGG)",
  "lac_vs_uhgg",   "Lactobacillus × Gut (UHGG)"
)

# Clamp extreme benefits so the hexbin stays readable
CLAMP <- 1.0

load_one <- function(sys_dir, label) {
  rows <- list()
  for (l in seq_len(nrow(SEL))) {
    p <- file.path(RESULTS, sys_dir, paste0(SEL$level_file[l], ".tsv"))
    if (!file.exists(p)) next
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    v <- df |>
      filter(growth_a_alone > THRESHOLD, growth_b_alone > THRESHOLD) |>
      transmute(system = label,
                short = SEL$short[l],
                ba = pmin(pmax(benefit_a, -CLAMP), CLAMP),
                bb = pmin(pmax(benefit_b, -CLAMP), CLAMP))
    rows[[length(rows)+1]] <- v
  }
  bind_rows(rows)
}

dat <- bind_rows(map2(SYSTEMS$sys_dir, SYSTEMS$label, load_one)) |>
  mutate(system = factor(system, levels = SYSTEMS$label),
         short  = factor(short, levels = SEL$short))

theme_gm <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title       = element_text(face = "bold", size = 7.5, hjust = 0),
      strip.text       = element_text(face = "bold", size = 6),
      strip.background = element_blank(),
      axis.title       = element_text(size = 6.5),
      axis.text        = element_text(size = 5.5, colour = "black"),
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25, colour = "black"),
      legend.position  = "right",
      legend.title     = element_text(size = 6),
      legend.text      = element_text(size = 5.5),
      legend.key.width = unit(0.25, "cm"),
      legend.key.height= unit(0.5, "cm"),
      plot.margin      = margin(2, 4, 2, 2)
    )
}

quad_lines <- list(
  geom_hline(yintercept = c(-EPS_BETA, EPS_BETA), colour = "grey55",
             linewidth = 0.2, linetype = "dashed"),
  geom_vline(xintercept = c(-EPS_BETA, EPS_BETA), colour = "grey55",
             linewidth = 0.2, linetype = "dashed")
)

figS3 <- ggplot(dat, aes(ba, bb)) +
  geom_bin2d(bins = 40) +
  quad_lines +
  facet_grid(system ~ short) +
  scale_fill_viridis_c(trans = "log10", option = "magma",
                       name = "Pairs", labels = comma) +
  scale_x_continuous(limits = c(-CLAMP, CLAMP), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-CLAMP, CLAMP), breaks = c(-1, 0, 1)) +
  labs(title = expression(paste(bold("Benefit landscape "),
                                 beta[A], bold(" vs "), beta[B],
                                 bold(" across the gradient"))),
       x = expression(beta[A]~"(probiotic relative benefit)"),
       y = expression(beta[B]~"(commensal relative benefit)")) +
  theme_gm()

ggsave(file.path(OUTDIR, "figS3_benefit_landscape.tiff"), figS3,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "figS3_benefit_landscape.pdf"), figS3,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "figS3_benefit_landscape.png"), figS3,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved figS3_benefit_landscape (TIFF/PDF/PNG) to", OUTDIR, "\n")
