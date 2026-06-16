#!/usr/bin/env Rscript
# =============================================================================
# plot_fig2_phylo_trees.R   (Main Figure 2, phylogeny panels)
# -----------------------------------------------------------------------------
# Visualise the de novo bac120 trees produced by GTDB-Tk (de_novo_wf) / IQ-TREE
# for the selected probiotic strains. One panel per genus, coloured by species.
#
# Inputs (newick, produced by run_gtdbtk_denovo.sh):
#   test/akk/iqtree_out/akk_bac120.treefile   (or the GTDB-Tk *decorated.tree)
#   test/lac/iqtree_out/lac_bac120.treefile
#
# Install once:
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("ggtree", "treeio"))
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtree)     # Bioconductor
  library(treeio)     # Bioconductor
  library(patchwork)
})

AKK_TREE <- "test/akk/iqtree_out/akk_bac120.treefile"
LAC_TREE <- "test/lac/iqtree_out/lac_bac120.treefile"
OUTDIR   <- "results/figures_paper"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── tip-label → pretty name + species group (edit to match your accessions) ──
akk_meta <- tribble(
  ~tip,                     ~label,                        ~species,
  "akk_GCF_000020225.1",    "A. muciniphila ATCC BAA-835", "A. muciniphila",
  "akk_GCF_018847155.1",    "A. muciniphila II",           "A. muciniphila",
  "akk_GCF_018847135.1",    "A. massiliensis I",           "A. massiliensis",
  "akk_GCF_023516715.1",    "A. massiliensis II",          "A. massiliensis",
  "akk_GCF_018847015.1",    "A. biwaensis I",              "A. biwaensis",
  "akk_GCF_026072915.1",    "A. biwaensis II",             "A. biwaensis",
  "akk_GCF_001683795.1",    "A. glycaniphila (outgroup)",  "outgroup"
)

lac_meta <- tribble(
  ~tip,                                        ~label,                 ~species,
  "L_acidophilus_NCFM_GCF_000011985.1",        "L. acidophilus NCFM",  "Lactobacillus",
  "L_crispatus_ST1_GCF_000091765.1",           "L. crispatus ST1",     "Lactobacillus",
  "L_crispatus_M247_GCF_026740115.1",          "L. crispatus M247",    "Lactobacillus",
  "L_gasseri_ATCC33323_GCF_000014425.1",       "L. gasseri ATCC33323", "Lactobacillus",
  "L_delbrueckii_ATCC_BAA365_GCF_000014405.1", "L. delbrueckii",       "Lactobacillus",
  "L_casei_ATCC393_GCF_000829055.1",           "L. casei ATCC393",     "Lacticaseibacillus",
  "L_paracasei_ATCC334_GCF_000014525.1",       "L. paracasei ATCC334", "Lacticaseibacillus",
  "L_rhamnosus_GG_GCF_000026505.1",            "L. rhamnosus GG",      "Lacticaseibacillus",
  "L_plantarum_WCFS1_GCF_000203855.3",         "L. plantarum WCFS1",   "Lactiplantibacillus",
  "L_reuteri_DSM20016_GCF_000016825.1",        "L. reuteri DSM20016",  "Limosilactobacillus",
  "lac_GCF_000014445.1",                       "Leuconostoc (outgroup)","outgroup"
)

SPECIES_COLOURS <- c(
  "A. muciniphila"      = "#762a83",
  "A. massiliensis"     = "#9970ab",
  "A. biwaensis"        = "#c2a5cf",
  "Lactobacillus"       = "#1b7837",
  "Lacticaseibacillus"  = "#5aae61",
  "Lactiplantibacillus" = "#a6dba0",
  "Limosilactobacillus" = "#d9f0d3",
  "outgroup"            = "grey60"
)

# ── helper: read tree, join metadata, draw ───────────────────────────────────
draw_tree <- function(tree_path, meta, title) {
  if (!file.exists(tree_path)) {
    warning("Tree not found: ", tree_path, " — run run_gtdbtk_denovo.sh first.")
    return(NULL)
  }
  tr <- read.tree(tree_path)
  # Build the base ggtree, then attach metadata WITHOUT clobbering the tree's
  # own `label` column (which holds tip names AND internal-node support values).
  # We add pretty names / species under distinct column names via match().
  p <- ggtree(tr, size = 0.35)
  idx <- match(p$data$label, meta$tip)
  p$data$pretty  <- meta$label[idx]
  p$data$species <- meta$species[idx]

  p +
    geom_tippoint(aes(colour = species), size = 1.6, na.rm = TRUE) +
    geom_tiplab(aes(label = pretty, colour = species),
                size = 2.0, family = "Arial", hjust = -0.02, na.rm = TRUE) +
    # bootstrap / support values on internal nodes (IQ-TREE stores them in label)
    geom_nodelab(aes(label = label), size = 1.4, hjust = 1.2, vjust = -0.4,
                 colour = "grey40", na.rm = TRUE) +
    scale_colour_manual(values = SPECIES_COLOURS, name = NULL, na.value = "black") +
    labs(title = title) +
    theme_tree() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 7.5, family = "Arial"),
          plot.margin = margin(2, 30, 2, 2)) +    # right margin for long tip labels
    coord_cartesian(clip = "off") +
    ggtree::geom_treescale(width = 0.1, fontsize = 1.6, linesize = 0.3,
                           offset = 0.3)
}

pAkk <- draw_tree(AKK_TREE, akk_meta, "Akkermansia (6 strains, 3 species)")
pLac <- draw_tree(LAC_TREE, lac_meta, "Lactobacillus group (10 strains, 5 genera)")

# Assemble as one column; combine with other Fig 1 panels in your layout
fig2_phylo <- (pAkk / pLac) + plot_layout(heights = c(1, 1.4))

ggsave(file.path(OUTDIR, "fig2_phylo_trees.tiff"), fig2_phylo,
       width = 90/25.4, height = 150/25.4, dpi = 300, device = "tiff", compression = "lzw")
ggsave(file.path(OUTDIR, "fig2_phylo_trees.pdf"), fig2_phylo,
       width = 90/25.4, height = 150/25.4, device = cairo_pdf)
ggsave(file.path(OUTDIR, "fig2_phylo_trees.png"), fig2_phylo,
       width = 90/25.4, height = 150/25.4, dpi = 300)

cat("Saved fig2_phylo_trees (TIFF/PDF/PNG) to", OUTDIR, "\n")
