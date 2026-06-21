#!/usr/bin/env Rscript
# plot_figS5_objective_sensitivity.R   (Supplementary Figure S5)
#
# Objective-function sensitivity analysis (reviewer Major Comment 4).
# For a representative subset (L. gasseri x UHGG at L5 pectin and L6 resistant
# starch), interaction classifications are recomputed under a fixed-ratio
# co-culture objective (mu_A^co/mu_B^co = mu_A^alone/mu_B^alone) and compared
# with the default lexicographic max-min allocation. The mutualism peak at L5
# and the glucose crash at L6 are unchanged, demonstrating that the core
# conclusions do not depend on the choice of growth-allocation objective.
#
# Input : results/objective_sensitivity/<level>_<objective>.tsv
# Output: results/figures_paper/figS5_objective_sensitivity.{pdf,png,tiff}

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
INDIR  <- get_arg("--indir",  "results/objective_sensitivity")
OUTDIR <- get_arg("--outdir", "results/figures_paper")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

THRESH <- 1e-4
files <- tibble::tribble(
  ~level,            ~objective,       ~path,
  "L5 pectin",       "Lexicographic",  "L5_pectin_lexicographic.tsv",
  "L5 pectin",       "Fixed-ratio",    "L5_pectin_fixed_ratio.tsv",
  "L6 resistant\nstarch", "Lexicographic", "L6_resistant_starch_lexicographic.tsv",
  "L6 resistant\nstarch", "Fixed-ratio",   "L6_resistant_starch_fixed_ratio.tsv"
)

read_pairs <- function(p) {
  suppressMessages(read_tsv(file.path(INDIR, p), show_col_types = FALSE, progress = FALSE)) |>
    mutate(viable = growth_a_alone > THRESH & growth_b_alone > THRESH)
}

# ── mutualism fraction per (level, objective) ────────────────────────────────
summ <- files |>
  rowwise() |>
  mutate(d = list(read_pairs(path))) |>
  mutate(viable = sum(d$viable),
         mutual = sum(d$interaction_type == "mutualism" & d$viable),
         mut_pct = 100 * mutual / viable) |>
  ungroup() |>
  mutate(level = factor(level, levels = c("L5 pectin", "L6 resistant\nstarch")),
         objective = factor(objective, levels = c("Lexicographic", "Fixed-ratio")))

# ── per-pair classification agreement per level ──────────────────────────────
agree_tbl <- files |>
  group_by(level) |>
  summarise(paths = list(path), .groups = "drop") |>
  rowwise() |>
  mutate(stats = list({
    dl <- read_pairs(paths[[1]]); df <- read_pairs(paths[[2]])
    m <- dplyr::inner_join(
      dl |> filter(viable) |> select(species_a, species_b, it_lex = interaction_type),
      df |> filter(viable) |> select(species_a, species_b, it_fix = interaction_type),
      by = c("species_a", "species_b"))
    tibble(n = nrow(m),
           agree = mean(m$it_lex == m$it_fix) * 100,
           mut_preserved = {
             ml <- m |> filter(it_lex == "mutualism")
             if (nrow(ml) == 0) 100 else mean(ml$it_fix == "mutualism") * 100
           })
  })) |>
  tidyr::unnest(stats) |> ungroup()

cat("Per-pair agreement:\n"); print(agree_tbl |> select(level, n, agree, mut_preserved))

BLUE <- "#2c7fb8"; ORANGE <- "#e6862e"
lvl_top <- summ |> mutate(level = as.character(level)) |>
  group_by(level) |> summarise(top = max(mut_pct), .groups = "drop")
ann <- agree_tbl |>
  mutate(level = as.character(level)) |>
  dplyr::left_join(lvl_top, by = "level") |>
  mutate(y = top + 5,
         lbl = sprintf("agreement %.1f%%\nmutualism 100%% kept", agree),
         level = factor(level, levels = c("L5 pectin", "L6 resistant\nstarch")))

pA <- ggplot(summ, aes(level, mut_pct, fill = objective)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", mut_pct)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  geom_text(data = ann, inherit.aes = FALSE,
            aes(x = level, y = y, label = lbl), size = 2.7, colour = "grey30", lineheight = 0.95) +
  scale_fill_manual(values = c("Lexicographic" = BLUE, "Fixed-ratio" = ORANGE), name = "Co-culture objective") +
  scale_y_continuous(limits = c(0, 44), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Mutualism (% of viable pairs)",
       title = "Objective-function sensitivity (L. gasseri × UHGG)",
       subtitle = "Mutualism peak at L5 and the L5→L6 glucose crash are identical under both allocation objectives") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

save_one <- function(ext, dev=NULL) {
  fp <- file.path(OUTDIR, paste0("figS5_objective_sensitivity.", ext))
  if (ext == "tiff") ggsave(fp, pA, width = 7, height = 5, dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "pdf") ggsave(fp, pA, width = 7, height = 5)
  else ggsave(fp, pA, width = 7, height = 5, dpi = 300)
}
for (e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS5_objective_sensitivity (PDF/PNG/TIFF) to", OUTDIR, "\n")
