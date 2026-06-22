#!/usr/bin/env Rscript
# plot_figS7_glucose_isolation.R   (Supplementary Figure S7)
#
# Controlled carbon-quality contrast (reviewer Major Comment 2a).
# Lactobacillus × UHGG, all run with the same binary. The key contrast holds the
# L1–L4 carbon background IDENTICAL and changes only the last carbon source added:
#   L5      = L4 background + pectin   (structurally complex)   — mutualism peak
#   L5-glc  = L4 background + glucose  (easily fermentable sugar, SAME background)
# Replacing pectin with glucose, against an identical background, lowers mutualism
# (21.7 % -> 18.4 %): a clean test that carbon QUALITY — not richness — sets the
# direction of cooperation. L6 (which layers further free sugars: glucose, maltose,
# maltodextrin on top of L5) is shown for reference: more easily fermentable sugar
# deepens the crash to 16.6 %.
# (Akkermansia is omitted: lacking general sugar transporters, it does not grow on
#  a glucose-substituted medium.)
#
# Input : results/glucose_isolation/lac_{L5_pectin,L5glc,L6_resistant_starch}.tsv
# Output: results/figures_paper/figS7_glucose_isolation.{pdf,png,tiff}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
INDIR  <- get_arg("--indir",  "results/glucose_isolation")
OUTDIR <- get_arg("--outdir", "results/figures_paper")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
THRESH <- 1e-4

spec <- tibble::tribble(
  ~level,                          ~file,                         ~grp,
  "L5\nL4 + pectin",               "lac_L5_pectin",               "complex carbon",
  "L5-glc\nL4 + glucose",          "lac_L5glc",                   "free glucose",
  "L6\nL5 + more free sugars",     "lac_L6_resistant_starch",     "more free sugars"
)

mut_pct <- function(file) {
  d <- suppressMessages(read_tsv(file.path(INDIR, paste0(file, ".tsv")),
                                 show_col_types = FALSE, progress = FALSE)) |>
    mutate(viable = growth_a_alone > THRESH & growth_b_alone > THRESH)
  tibble(viable = sum(d$viable),
         mutual = sum(d$interaction_type == "mutualism" & d$viable),
         mut_pct = 100 * sum(d$interaction_type == "mutualism" & d$viable) / sum(d$viable))
}

summ <- spec |> rowwise() |> mutate(s = list(mut_pct(file))) |>
  tidyr::unnest(s) |> ungroup() |>
  mutate(level = factor(level, levels = spec$level),
         grp   = factor(grp, levels = c("complex carbon", "free glucose", "more free sugars")))

cat("Lactobacillus × UHGG mutualism (same binary):\n")
print(summ |> select(level, viable, mutual, mut_pct))

cols <- c("complex carbon" = "#41ab5d", "free glucose" = "#e6862e", "more free sugars" = "#d6604d")

p <- ggplot(summ, aes(level, mut_pct, fill = grp)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", mut_pct)), vjust = -0.6, size = 4, fontface = "bold") +
  # bracket marking the controlled contrast L5 -> L5-glc
  annotate("segment", x = 1, xend = 2, y = 23.4, yend = 23.4, colour = "grey35", linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = 23.0, yend = 23.4, colour = "grey35", linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = 23.0, yend = 23.4, colour = "grey35", linewidth = 0.4) +
  annotate("text", x = 1.5, y = 24.1,
           label = "same L1–L4 background\npectin → glucose:  −3.3 points",
           size = 2.9, colour = "grey25", lineheight = 0.9) +
  scale_fill_manual(values = cols, name = "Last carbon added") +
  scale_y_continuous(limits = c(0, 26), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Mutualism (% of viable pairs)",
       title = "Carbon quality sets the direction: glucose suppresses cooperation more than pectin",
       subtitle = "Same L1–L4 background: replacing pectin (L5) with glucose (L5-glc) lowers mutualism; further free sugars (L6) deepen the crash") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 9, lineheight = 0.9))

save_one <- function(ext) {
  fp <- file.path(OUTDIR, paste0("figS7_glucose_isolation.", ext))
  if (ext == "tiff") ggsave(fp, p, width = 7.6, height = 5, dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "pdf") ggsave(fp, p, width = 7.6, height = 5)
  else ggsave(fp, p, width = 7.6, height = 5, dpi = 300)
}
for (e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS7_glucose_isolation (PDF/PNG/TIFF) to", OUTDIR, "\n")
