#!/usr/bin/env Rscript
# plot_figS11_gpr_coverage.R  (Supplementary Figure S11)
# Cross-feeding gene support reflects transporter annotation, not a gap-filling artefact.
#
# Input : results/figS11/gpr_coverage.tsv
# Output: results/figures_paper/figS11_gpr_coverage.{pdf,png,tiff}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr) })

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
INFILE <- get_arg("--in", "results/figS11/gpr_coverage.tsv")
OUTDIR <- get_arg("--outdir", "results/figures_paper")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <-  80 / 25.4


d <- suppressMessages(read_tsv(INFILE, show_col_types = FALSE)) |>
  mutate(category = factor(category,
           levels = c("Intracellular reactions","Transport reactions","Cross-feeding flux")),
         panel = factor(panel, levels = c("Akkermansia","Lactobacillus")))

cols <- c("Intracellular reactions" = "#6baed6",
          "Transport reactions"     = "#fdae6b",
          "Cross-feeding flux"      = "#74c476")

p <- ggplot(d, aes(panel, gene_support_pct, fill = category)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", gene_support_pct)),
            position = position_dodge(width = 0.78), vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_x_discrete(labels = sp_md) +
  scale_y_continuous(limits = c(0, 92), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Gene-supported reactions / flux (%)",
       title = "Cross-feeding gene support reflects transporter annotation, not a gap-filling artefact",
       subtitle = "Intracellular reactions ~80% gene-supported; transport reactions lower; cross-feeding matches the transport level") +
  theme_classic(base_size = 8, base_family = BMC_FONT) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 8.5),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        legend.text = element_text(size = 6.5), legend.key.size = unit(0.32, "cm"),
        plot.subtitle = element_text(size = 6, colour = "grey30"),
        axis.text.x = md(size = 7))

save_one <- function(ext) {
  fp <- file.path(OUTDIR, paste0("figS11_gpr_coverage.", ext))
  if (ext == "tiff") ggsave(fp, p, width = FIG_W, height = FIG_H, dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "pdf") ggsave(fp, p, width = FIG_W, height = FIG_H, device = cairo_pdf)
  else ggsave(fp, p, width = FIG_W, height = FIG_H, dpi = 300)
}
for (e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS11_gpr_coverage (PDF/PNG/TIFF) to", OUTDIR, "\n")
