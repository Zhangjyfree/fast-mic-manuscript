#!/usr/bin/env Rscript
# =============================================================================
# plot_figS3_crossfeed_landscape.R   (Supplementary Figure S3)
# -----------------------------------------------------------------------------
# Distinct metabolic currencies of cooperation. For each system, the most
# frequently exchanged metabolites among mutualistic pairs (fraction of
# cooperating pairs that trade each compound). Supports the main-text claim
# that the generalist (Lactobacillus) and specialist (Akkermansia) cooperate
# through systematically different metabolite sets.
#
# Input: results/figures_paper/crossfeed_landscape_table.tsv
#        columns: sys_label, metabolite, cpd, n_pairs, n_mut, frac
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
TOP_N   <- as.integer(get_arg("--top", "15"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FIG_W <- 170 / 25.4
FIG_H <- 95 / 25.4
FIG_DPI <- 300

tbl_path <- file.path(OUTDIR, "crossfeed_landscape_table.tsv")
if (!file.exists(tbl_path)) stop("Missing input: ", tbl_path)

SYS_LEVELS <- c("Akkermansia × UHGG (gut)",
                "Lactobacillus × UHGG (gut)")
SYS_COLOURS <- setNames(c("#762a83", "#1b7837"), SYS_LEVELS)

cf <- suppressMessages(read_tsv(tbl_path, show_col_types = FALSE, progress = FALSE)) |>
  filter(sys_label %in% SYS_LEVELS) |>
  mutate(sys_label = factor(sys_label, levels = SYS_LEVELS))

# Top-N metabolites per system, ordered within facet by fraction.
# Manual within-facet ordering (avoids the tidytext dependency): make a unique
# per-facet label and set its factor levels by ascending frac within system.
top_cf <- cf |>
  group_by(sys_label) |>
  slice_max(frac, n = TOP_N, with_ties = FALSE) |>
  ungroup() |>
  arrange(sys_label, frac) |>
  mutate(met_key = paste(as.integer(sys_label), metabolite, sep = "___"),
         met_key = factor(met_key, levels = met_key))

theme_gm <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title       = element_text(face = "bold", size = 7.5, hjust = 0),
      strip.text       = element_text(face = "bold", size = 6.5),
      strip.background = element_blank(),
      axis.title       = element_text(size = 6.5),
      axis.text        = element_text(size = 5.5, colour = "black"),
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.15),
      legend.position  = "none",
      plot.margin      = margin(2, 4, 2, 2)
    )
}

figS3 <- ggplot(top_cf, aes(frac, met_key, fill = sys_label)) +
  geom_col(width = 0.78) +
  geom_text(aes(label = sprintf("%.0f%%", frac * 100)),
            hjust = -0.15, size = 1.7, family = "Arial", colour = "grey20") +
  facet_wrap(~ sys_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = SYS_COLOURS) +
  scale_y_discrete(labels = function(x) sub("^[0-9]+___", "", x)) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Fraction of mutualistic pairs exchanging metabolite",
       y = NULL,
       title = paste0("Top ", TOP_N,
                      " cross-fed metabolites in cooperating pairs, by system")) +
  theme_gm()

ggsave(file.path(OUTDIR, "figS3_crossfeed_landscape.tiff"), figS3,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "figS3_crossfeed_landscape.pdf"), figS3,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "figS3_crossfeed_landscape.png"), figS3,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved figS3_crossfeed_landscape (TIFF/PDF/PNG) to", OUTDIR, "\n")
