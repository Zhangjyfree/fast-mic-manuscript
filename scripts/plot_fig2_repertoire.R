#!/usr/bin/env Rscript
# Fig 2 v2: phylogeny + metabolic repertoire + niche breadth (tip-aligned).
# Per system: ML tree | trait bars (reactions/metabolites/genes/transporters) | growth heatmap (L0-L9).
# Operationalises "generalist (Lactobacillus) vs specialist (Akkermansia)" as data.
#
# Inputs (from scripts/extract_fig2_traits.py):
#   results/fig2/fig2_traits_{sys}.tsv , results/fig2/fig2_growth_{sys}.tsv
#   trees: test/{sys}/iqtree_out/{sys}_bac120.treefile
# Usage: Rscript scripts/plot_fig2_v2.R [output_prefix]

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(ggtree); library(treeio); library(aplot); library(patchwork); library(scales); library(cowplot)
})

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

options(lifecycle_verbosity = "quiet")  # silence ggtree's internal ggplot2 deprecation warnings

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[1] else "results/figures_paper/fig2_v2_repertoire"
FM <- "."   # repo root

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 210 / 25.4      # within the 225 mm limit

LEVELS <- c("L0","L1","L2","L3","L4","L5","L6","L7","L8","L9")
TRAITS <- c("reactions","transporters","exchanges")
TRAIT_LAB <- c(reactions="Reactions", transporters="Transporters", exchanges="Exchanges")

# akk species by accession (manuscript methods); lac genus by name prefix
akk_species <- function(s) {
  dplyr::case_when(
    grepl("000020225|018847155", s) ~ "A. muciniphila",
    grepl("018847135|023516715", s) ~ "A. massiliensis",
    grepl("018847015|026072915", s) ~ "A. biwaensis",
    TRUE ~ "Akkermansia sp.")
}
lac_genus <- function(s) {
  dplyr::case_when(
    grepl("casei|paracasei|rhamnosus", s, ignore.case = TRUE) ~ "Lacticaseibacillus",
    grepl("plantarum", s, ignore.case = TRUE)                 ~ "Lactiplantibacillus",
    grepl("reuteri", s, ignore.case = TRUE)                   ~ "Limosilactobacillus",
    TRUE ~ "Lactobacillus")
}
pretty_lac <- function(s) s |> sub("_GCF_.*$", "", x = _) |> gsub("_", " ", x = _) |>
  sub("^([A-Z]) ", "\\1. ", x = _)   # "L acidophilus" -> "L. acidophilus"
pretty_akk <- function(s) paste0(akk_species(s), " ", sub("^akk_", "", s))  # species + accession
grp_cols <- c("A. muciniphila"="#762A83","A. massiliensis"="#9970AB","A. biwaensis"="#C2A5CF",
              "Lactobacillus"="#1B7837","Lacticaseibacillus"="#5AAE61",
              "Lactiplantibacillus"="#A6DBA0","Limosilactobacillus"="#00695C")

base_theme <- theme_bw(base_size = 7, base_family = BMC_FONT) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 6.5),
        axis.text = element_text(size = 5.5), legend.position = "bottom",
        legend.key.size = unit(0.28, "cm"), legend.text = md(size = 6),
        legend.title = element_text(size = 6.5), strip.text = element_text(size = 6),
        strip.background = element_rect(fill = "grey92", colour = NA),
        plot.title = element_text(size = 8, face = "bold"))

build_system <- function(sys_key, group_fun, pretty_fun, title, growth_lim) {
  tree <- read.tree(sprintf("%s/test/%s/iqtree_out/%s_bac120.treefile", FM, sys_key, sys_key))
  traits <- read_tsv(sprintf("%s/results/fig2/fig2_traits_%s.tsv", FM, sys_key), show_col_types = FALSE)
  growth <- read_tsv(sprintf("%s/results/fig2/fig2_growth_%s.tsv", FM, sys_key), show_col_types = FALSE)

  # keep only strains present in the tree (drop outgroup tips not in trait table)
  keep <- intersect(tree$tip.label, traits$strain)
  tree <- ape::keep.tip(tree, keep)
  meta <- tibble(strain = tree$tip.label) |>
    mutate(group = group_fun(strain), pretty = pretty_fun(strain),
           pretty_it = sp_plotmath(pretty))

  # ── tree (carries the system title; leftmost, so no title collision) ──
  ptree <- ggtree(tree, linewidth = 0.4) %<+% meta +
    geom_tippoint(aes(colour = group), size = 1.2) +
    geom_tiplab(aes(label = pretty_it, colour = group), parse = TRUE,
                family = BMC_FONT, size = 1.85, align = TRUE,
                linetype = "dotted", linesize = 0.2, offset = 0.005, show.legend = FALSE) +
    scale_colour_manual(values = grp_cols, name = NULL, labels = sp_md) +
    ggtree::geom_treescale(width = 0.05, fontsize = 1.5, linesize = 0.25, offset = 0.4,
                           family = BMC_FONT) +
    ## aplot mis-measures a mixed italic/roman run, so the genus goes in the
    ## title (italic throughout) and its role in the subtitle (roman throughout)
    labs(title = sp_md(sub(" \\(.*$", "", title)),
         subtitle = sub("^.*\\((.*)\\)$", "\\1", title)) +
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2),  # colours shown via labelled tips
          plot.title = md(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.5, colour = "grey35",
                                       family = BMC_FONT)) +
    xlim(0, max(ape::node.depth.edgelength(tree)) * 3.4)   # headroom for aligned tip labels

  # ── trait bars (faceted, y = strain) ──
  tl <- traits |>
    select(strain, all_of(TRAITS)) |>
    pivot_longer(-strain, names_to = "trait", values_to = "value") |>
    mutate(trait = factor(TRAIT_LAB[trait], levels = TRAIT_LAB[TRAITS]),
           strain = factor(strain, levels = rev(get_taxa_name(ptree))))
  ptraits <- ggplot(tl, aes(value, strain, fill = trait)) +
    geom_col(width = 0.66, show.legend = FALSE) +
    geom_text(aes(label = value), hjust = -0.08, size = 1.7, colour = "grey25") +
    facet_wrap(~ trait, nrow = 1, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
    scale_fill_manual(values = c("Reactions"="#4393C3","Transporters"="#D6604D","Exchanges"="#5AAE61")) +
    labs(x = NULL, y = NULL, title = NULL) +
    base_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())  # values labelled on bars

  # ── growth heatmap (x = level, y = strain) ──
  gh <- growth |>
    mutate(level = factor(level, levels = LEVELS),
           strain = factor(strain, levels = rev(get_taxa_name(ptree))))
  pgrow <- ggplot(gh, aes(level, strain, fill = growth)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_viridis_c(option = "D", limits = growth_lim, guide = "none") +  # shared colourbar drawn separately
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = NULL) +
    base_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(size = 5.2))

  # align by tip: tree | traits | growth  (growth on the far right, next to its colourbar)
  ptraits |> insert_left(ptree, width = 0.85) |> insert_right(pgrow, width = 0.62)
}

# shared growth colour scale across both systems (comparable A vs B)
grng <- range(c(read_tsv(sprintf("%s/results/fig2/fig2_growth_akk.tsv", FM), show_col_types = FALSE)$growth,
                read_tsv(sprintf("%s/results/fig2/fig2_growth_lac.tsv", FM), show_col_types = FALSE)$growth),
              na.rm = TRUE)
# shared vertical growth colourbar, extracted from a source plot and placed on the far right
cbar <- cowplot::get_legend(
  ggplot(data.frame(x = 1, growth = grng), aes(x, x, fill = growth)) + geom_tile() +
    scale_fill_viridis_c(option = "D", limits = grng,
      name = "Monoculture\ngrowth (h\u207b\u00b9)",
      guide = guide_colourbar(barwidth = unit(0.3, "cm"), barheight = unit(3.2, "cm"))) +
    theme(text = element_text(family = BMC_FONT),
          legend.position = "right",
          legend.title = element_text(size = 6.5, family = BMC_FONT),
          legend.text = element_text(size = 5.5, family = BMC_FONT)))

suppressWarnings({   # ggtree internals emit benign check_dots/deprecation warnings at build & draw
  gakk <- build_system("akk", akk_species, pretty_akk, "Akkermansia (specialist)", grng)
  glac <- build_system("lac", lac_genus,   pretty_lac, "Lactobacillus-group (generalist)", grng)
  stacked  <- aplot::plot_list(gakk, glac, ncol = 1, heights = c(6, 10), tag_levels = "A") &
    theme(plot.tag = element_text(family = BMC_FONT, face = "bold"))
  combined <- cowplot::plot_grid(stacked, cbar, ncol = 2, rel_widths = c(1, 0.11))
  for (ext in c("pdf","png","tiff")) {
    path <- paste0(out_prefix, ".", ext)
    if (ext == "tiff") ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in",
                              dpi = 300, device = "tiff", compression = "lzw")
    else if (ext == "png") ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in", dpi = 300)
    else ggsave(path, combined, width = FIG_W, height = FIG_H, units = "in",
                device = cairo_pdf)   # cairo: keeps Arial + plotmath italics
    cat("Saved:", path, "\n")
  }
})
