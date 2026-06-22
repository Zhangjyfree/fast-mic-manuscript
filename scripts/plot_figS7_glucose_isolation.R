#!/usr/bin/env Rscript
# plot_figS7_glucose_isolation.R   (Supplementary Figure S7)
#
# Controlled carbon-quality contrast, BOTH probiotic systems (reviewer Major
# Comment 2a; directly answers Question (ii): is the cooperation collapse shared
# by the specialist AND the generalist?). The L1-L4 carbon background is held
# IDENTICAL and only the last carbon source is changed:
#   L5      = L4 background + pectin   (structurally complex)
#   L5-glc  = L4 background + glucose  (easily fermentable sugar, SAME background)
#   L6      = L5 + further free sugars (glucose, maltose, maltodextrin), reference
# Replacing pectin with glucose, against an identical background, lowers mutualism
# in BOTH systems:
#   Akkermansia   (specialist) 10.4 % -> 6.5 %  (-3.9)
#   Lactobacillus (generalist) 21.7 % -> 18.4 % (-3.3)
# i.e. the glucose-driven suppression of cooperation is a carbon-QUALITY effect
# (not richness) and it is SHARED by both probiotics. Note: Akkermansia itself
# cannot import glucose; the effect propagates through its UHGG partners, which
# can — easy sugar lets partners self-supply, eroding the cross-feeding niche.
#
# Input : results/glucose_isolation/{akk,lac}_{L5_pectin,L5glc,L6_resistant_starch}.tsv
# Output: results/figures_paper/figS7_glucose_isolation.{pdf,png,tiff}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
INDIR  <- get_arg("--indir",  "results/glucose_isolation")
OUTDIR <- get_arg("--outdir", "results/figures_paper")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
THRESH <- 1e-4

spec <- tidyr::expand_grid(
  tibble::tribble(~system,                      ~pre,
                  "Akkermansia (specialist)",   "akk",
                  "Lactobacillus (generalist)", "lac"),
  tibble::tribble(~level,                       ~file,                     ~grp,
                  "L5\nL4 + pectin",            "L5_pectin",               "complex carbon",
                  "L5-glc\nL4 + glucose",       "L5glc",                   "free glucose",
                  "L6\n+ more free sugars",     "L6_resistant_starch",     "more free sugars"))

mut_pct <- function(pre, file) {
  d <- suppressMessages(read_tsv(file.path(INDIR, paste0(pre, "_", file, ".tsv")),
                                 show_col_types = FALSE, progress = FALSE)) |>
    mutate(viable = growth_a_alone > THRESH & growth_b_alone > THRESH)
  100 * sum(d$interaction_type == "mutualism" & d$viable) / sum(d$viable)
}

summ <- spec |> rowwise() |> mutate(mut_pct = mut_pct(pre, file)) |> ungroup() |>
  mutate(system = factor(system, levels = c("Akkermansia (specialist)", "Lactobacillus (generalist)")),
         level  = factor(level,  levels = unique(spec$level)),
         grp    = factor(grp,    levels = c("complex carbon", "free glucose", "more free sugars")))

cat("Mutualism % by system × medium (same binary):\n")
print(summ |> select(system, level, mut_pct), n = 99)

cols <- c("complex carbon" = "#41ab5d", "free glucose" = "#e6862e", "more free sugars" = "#d6604d")

# per-facet annotation for the controlled L5 -> L5-glc drop
ann <- summ |> filter(grp != "more free sugars") |>
  group_by(system) |>
  summarise(ytop = max(mut_pct) * 1.16,
            delta = mut_pct[level == levels(level)[2]] - mut_pct[level == levels(level)[1]],
            .groups = "drop") |>
  mutate(lab = sprintf("pectin → glucose:  %+.1f", delta))

p <- ggplot(summ, aes(level, mut_pct, fill = grp)) +
  geom_col(width = 0.64, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", mut_pct)), vjust = -0.55, size = 3.4, fontface = "bold") +
  geom_text(data = ann, aes(x = 1.5, y = ytop, label = lab),
            inherit.aes = FALSE, size = 3, colour = "grey25") +
  facet_wrap(~system, scales = "free_y") +
  scale_fill_manual(values = cols, name = "Last carbon added") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Mutualism (% of viable pairs)",
       title = "Carbon quality sets the direction in both probiotics",
       subtitle = paste0("Same L1–L4 background: replacing pectin (L5) with glucose (L5-glc) lowers mutualism in BOTH the specialist (−3.9)\n",
                         "and the generalist (−3.3) — the glucose-driven collapse of cooperation is shared, not system-specific")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text = element_text(face = "bold.italic", size = 9.5),
        axis.text.x = element_text(size = 8, lineheight = 0.9))

save_one <- function(ext) {
  fp <- file.path(OUTDIR, paste0("figS7_glucose_isolation.", ext))
  if (ext == "tiff") ggsave(fp, p, width = 8.6, height = 5, dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "pdf") ggsave(fp, p, width = 8.6, height = 5)
  else ggsave(fp, p, width = 8.6, height = 5, dpi = 300)
}
for (e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS7_glucose_isolation (PDF/PNG/TIFF) to", OUTDIR, "\n")
