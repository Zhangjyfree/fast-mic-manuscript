#!/usr/bin/env Rscript
# plot_figS8_gpr_coverage.R   (Supplementary Figure S8)
#
# Gap-filling control (reviewer Major Comment 3).
# Genome-wide GPR (gene-protein-reaction) coverage for the two probiotic panels,
# split into intracellular vs transport reactions, compared with the gene-supported
# fraction of cross-feeding flux. Intracellular reactions are well annotated (~80%),
# whereas transport reactions are systematically less annotated in all draft GEMs.
# The cross-feeding gene-support tracks the TRANSPORT level — not the genome-wide
# average — showing that the lower gene support of cross-feeding reflects the known
# difficulty of annotating transporters, NOT a gap-filling artefact that
# specifically inflates cooperation.
#
# Input : results/gpr_analysis/gpr_coverage.tsv
# Output: results/figures_paper/figS8_gpr_coverage.{pdf,png,tiff}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
INFILE <- get_arg("--in", "results/gpr_analysis/gpr_coverage.tsv")
OUTDIR <- get_arg("--outdir", "results/figures_paper")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

d <- suppressMessages(read_tsv(INFILE, show_col_types = FALSE)) |>
  mutate(category = factor(category,
           levels = c("Intracellular reactions","Transport reactions","Cross-feeding flux")),
         panel = factor(panel, levels = c("Akkermansia","Lactobacillus")))

cat("GPR coverage / cross-feeding gene-support:\n"); print(d)

cols <- c("Intracellular reactions" = "#6baed6",
          "Transport reactions"     = "#fdae6b",
          "Cross-feeding flux"      = "#74c476")

p <- ggplot(d, aes(panel, gene_support_pct, fill = category)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", gene_support_pct)),
            position = position_dodge(width = 0.78), vjust = -0.4, size = 3) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0, 92), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Gene-supported reactions / flux (%)",
       title = "Cross-feeding gene support reflects transporter annotation, not a gap-filling artefact",
       subtitle = "Intracellular reactions ~80% gene-supported; transport reactions lower; cross-feeding matches the transport level") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 10, face = "italic"))

save_one <- function(ext) {
  fp <- file.path(OUTDIR, paste0("figS8_gpr_coverage.", ext))
  if (ext == "tiff") ggsave(fp, p, width = 7.6, height = 5, dpi = 300, device = "tiff", compression = "lzw")
  else if (ext == "pdf") ggsave(fp, p, width = 7.6, height = 5)
  else ggsave(fp, p, width = 7.6, height = 5, dpi = 300)
}
for (e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS8_gpr_coverage (PDF/PNG/TIFF) to", OUTDIR, "\n")
