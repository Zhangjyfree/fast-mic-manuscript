#!/usr/bin/env Rscript
# plot_figS1_l0_substrates.R   (Supplementary Figure S1)
#
# L0 viability of the generalist panel (reviewer Major Comment 1).
# All ten Lactobacillus-group strains are viable in L0 monoculture. This figure
# shows, per strain, the uptake flux at L0 grouped by substrate class, revealing
# that growth is driven mainly by amino-acid fermentation, supplemented by
# mucin-derived amino sugars (N-acetylglucosamine, glucosamine, fucose, etc.);
# no strain uses fatty-acid oxidation. Uptake fluxes are from a parsimonious FBA
# under the fast-mic L0 medium; the monoculture growth rate annotation is the
# fast-mic value (h^-1).
#
# Input : results/figS1/lac_L0_uptake.tsv
# Output: results/figures_paper/figS1_l0_substrates.{pdf,png,tiff}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(tidyr) })

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
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
INFILE <- get_arg("--in", "results/figS1/lac_L0_uptake.tsv")
OUTDIR <- get_arg("--outdir", "results/figures_paper")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 105 / 25.4


# fast-mic L0 monoculture growth (h^-1) per strain
fm_growth <- c(
  "L. rhamnosus"=0.80, "L. casei"=0.79, "L. paracasei"=0.74, "L. plantarum"=0.73,
  "L. acidophilus"=0.42, "L. crispatus M247"=0.39, "L. gasseri"=0.39,
  "L. crispatus ST1"=0.39, "L. reuteri"=0.26, "L. delbrueckii"=0.14)

pretty <- function(s) {
  s <- sub("^L_", "L. ", s)
  s <- gsub("_GCF.*$", "", s)
  s <- sub("_NCFM|_ATCC393|_ATCC334|_WCFS1|_DSM20016|_GG|_ATCC_BAA365|_ATCC33323", "", s)
  s <- sub("_M247", " M247", s); s <- sub("_ST1", " ST1", s)
  s
}

d <- suppressMessages(read_tsv(INFILE, show_col_types = FALSE)) |>
  mutate(strain = pretty(strain),
         category = factor(category,
           levels = c("Amino acids","Mucin amino sugars","Other carbon")))
# order strains by fast-mic growth (desc)
ord <- names(sort(fm_growth, decreasing = TRUE))
d <- d |> mutate(strain = factor(strain, levels = ord))
lab <- tibble(strain = factor(ord, levels = ord), g = fm_growth[ord],
              tot = sapply(ord, function(s) sum(d$uptake_flux[d$strain==s])))

cols <- c("Amino acids"="#41ab5d", "Mucin amino sugars"="#e6862e", "Other carbon"="#bdbdbd")

p <- ggplot(d, aes(strain, uptake_flux, fill = category)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.2) +
  geom_text(data = lab, inherit.aes = FALSE,
            aes(strain, tot, label = sprintf("μ=%.2f", g)),
            vjust = -0.4, size = 1.9, colour = "grey25") +
  scale_fill_manual(values = cols, name = "Substrate class") +
  scale_x_discrete(labels = sp_md) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "L0 uptake flux (mmol gDW⁻¹ h⁻¹)",
       title = sp_md("At L0, Lactobacillus grows mainly by amino-acid fermentation"),
       subtitle = "All ten strains are viable in L0 (fast-mic μ = 0.14–0.80 h⁻¹); growth draws on amino acids and mucin amino sugars, not fatty-acid oxidation") +
  theme_classic(base_size = 8, base_family = BMC_FONT) +
  theme(legend.position = "top",
        plot.title = md(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 5.8, colour = "grey30"),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"),
        ## species italic, strain designation roman -> markdown, not face = "italic"
        axis.text.x = md(angle = 35, hjust = 1, size = 6.2))

save_one <- function(ext) {
  fp <- file.path(OUTDIR, paste0("figS1_l0_substrates.", ext))
  if (ext == "tiff") ggsave(fp, p, width = FIG_W, height = FIG_H, dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "pdf") ggsave(fp, p, width = FIG_W, height = FIG_H, device=cairo_pdf)
  else ggsave(fp, p, width = FIG_W, height = FIG_H, dpi = 300)
}
for (e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS1_l0_substrates (PDF/PNG/TIFF) to", OUTDIR, "\n")
