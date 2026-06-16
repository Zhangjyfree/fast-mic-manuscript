#!/usr/bin/env Rscript
# =============================================================================
# plot_fig6_strain_heatmap.R   (Main Figure 6)
# -----------------------------------------------------------------------------
# Strain-level heterogeneity in cooperation across the prebiotic gradient.
#
# Layout:
#   A  Heatmap: 10 Lactobacillus strains × 10 levels, faceted by partner pool
#   B  Per-strain trajectories (gut pool)
#   C  Akkermansia heatmap (6 strains, gut pool)
#
# Formatted for Gut Microbes (Taylor & Francis):
#   - Double-column: 170 mm, 300 DPI TIFF, Arial ≥6 pt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f,d){ i<-match(f,args); if(!is.na(i)&&i<length(args)) args[i+1] else d }
RESULTS <- get_arg("--results", "results")
OUTDIR  <- get_arg("--outdir",  file.path(RESULTS, "figures_paper"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

THRESHOLD <- 1e-4
FIG_W <- 170 / 25.4   # 6.69 in
FIG_H <- 155 / 25.4   # two stacked heatmaps (Lac, Akk)
FIG_DPI <- 300

LEVELS <- tibble::tribble(
  ~level_idx, ~level_file,              ~short,
  0, "L0_base",             "L0",
  1, "L1_inulin",           "L1",
  2, "L2_fos",              "L2",
  3, "L3_gos",              "L3",
  4, "L4_xos",              "L4",
  5, "L5_pectin",           "L5",
  6, "L6_resistant_starch", "L6",
  7, "L7_bglucan",          "L7",
  8, "L8_hmo",              "L8",
  9, "L9_mos",              "L9"
)

LAC_STRAINS <- tibble::tribble(
  ~raw,                                          ~display,                 ~group,
  "L_rhamnosus_GG_GCF_000026505.1",              "L. rhamnosus GG",       "Gut",
  "L_acidophilus_NCFM_GCF_000011985.1",          "L. acidophilus NCFM",   "Gut",
  "L_plantarum_WCFS1_GCF_000203855.3",           "L. plantarum WCFS1",    "Gut",
  "L_reuteri_DSM20016_GCF_000016825.1",          "L. reuteri DSM20016",   "Gut",
  "L_crispatus_ST1_GCF_000091765.1",             "L. crispatus ST1",      "Vaginal",
  "L_crispatus_M247_GCF_026740115.1",            "L. crispatus M247",     "Vaginal",
  "L_gasseri_ATCC33323_GCF_000014425.1",         "L. gasseri ATCC33323",  "Vaginal",
  "L_casei_ATCC393_GCF_000829055.1",             "L. casei ATCC393",      "Multi-site",
  "L_paracasei_ATCC334_GCF_000014525.1",         "L. paracasei ATCC334",  "Multi-site",
  "L_delbrueckii_ATCC_BAA365_GCF_000014405.1",   "L. delbrueckii",        "Dairy"
)

AKK_STRAINS <- tibble::tribble(
  ~raw,                         ~display,                       ~group,
  "akk_GCF_000020225.1",       "A. muciniphila ATCC BAA-835",  "muciniphila",
  "akk_GCF_018847155.1",       "A. muciniphila II",            "muciniphila",
  "akk_GCF_018847135.1",       "A. massiliensis I",            "massiliensis",
  "akk_GCF_023516715.1",       "A. massiliensis II",           "massiliensis",
  "akk_GCF_018847015.1",       "A. biwaensis I",               "biwaensis",
  "akk_GCF_026072915.1",       "A. biwaensis II",              "biwaensis"
)

GROUP_COLOURS <- c("Gut" = "#1b7837", "Vaginal" = "#e08214",
                   "Multi-site" = "#3182bd", "Dairy" = "#b35806")

# ── Compute per-strain mutualism fraction ────────────────────────────────────
compute_strain_frac <- function(sys_dir, strains, pool_label) {
  rows <- list()
  for (l in seq_len(nrow(LEVELS))) {
    p <- file.path(RESULTS, sys_dir, paste0(LEVELS$level_file[l], ".tsv"))
    if (!file.exists(p)) next
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    agg <- df |>
      filter(growth_a_alone > THRESHOLD, growth_b_alone > THRESHOLD) |>
      mutate(is_mut = tolower(interaction_type) == "mutualism") |>
      group_by(species_a) |>
      summarise(n_viable = n(), n_mut = sum(is_mut), .groups = "drop") |>
      mutate(frac = if_else(n_viable > 0, n_mut / n_viable, NA_real_),
             pool_label = pool_label,
             level_idx = LEVELS$level_idx[l])
    rows[[length(rows)+1]] <- agg
  }
  bind_rows(rows) |> inner_join(strains, by = c("species_a" = "raw"))
}

lac_all <- compute_strain_frac("lac_vs_uhgg", LAC_STRAINS, "Gut (UHGG)")

# Reorder strains by mean mutualism (highest on top)
lac_order <- lac_all |>
  group_by(display) |>
  summarise(mean_mut = mean(frac, na.rm = TRUE), .groups = "drop") |>
  arrange(mean_mut) |> pull(display)   # ascending → rev in factor = highest on top

lac_grid <- expand_grid(display = LAC_STRAINS$display,
                        pool_label = "Gut (UHGG)",
                        level_idx = LEVELS$level_idx) |>
  left_join(distinct(LAC_STRAINS, display, group), by = "display") |>
  left_join(lac_all |> select(display, pool_label, level_idx, n_viable, n_mut, frac),
            by = c("display","pool_label","level_idx")) |>
  left_join(LEVELS, by = "level_idx") |>
  mutate(display = factor(display, levels = lac_order),
         short = factor(short, levels = LEVELS$short),
         pool_label = factor(pool_label, levels = "Gut (UHGG)"),
         no_growth = is.na(frac))

akk_data <- compute_strain_frac("akk_vs_uhgg", AKK_STRAINS, "Gut (UHGG)")
akk_order <- akk_data |>
  group_by(display) |>
  summarise(mean_mut = mean(frac, na.rm = TRUE), .groups = "drop") |>
  arrange(mean_mut) |> pull(display)

akk_grid <- expand_grid(display = AKK_STRAINS$display,
                        level_idx = LEVELS$level_idx) |>
  left_join(distinct(AKK_STRAINS, display, group), by = "display") |>
  left_join(akk_data |> select(display, level_idx, n_viable, n_mut, frac),
            by = c("display","level_idx")) |>
  left_join(LEVELS, by = "level_idx") |>
  mutate(display = factor(display, levels = akk_order),
         short = factor(short, levels = LEVELS$short),
         no_growth = is.na(frac))

# ── Gut Microbes theme ───────────────────────────────────────────────────────
theme_hm <- function() {
  theme_minimal(base_size = 7, base_family = "Arial") +
    theme(
      plot.title       = element_text(face = "bold", size = 7.5, hjust = 0),
      plot.subtitle    = element_text(colour = "grey40", size = 6),
      axis.text.x      = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y      = element_text(size = 6),
      strip.text       = element_text(face = "bold", size = 7),
      panel.grid       = element_blank(),
      legend.title     = element_text(size = 6),
      legend.text      = element_text(size = 5.5),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width  = unit(0.25, "cm"),
      plot.tag         = element_text(face = "bold", size = 9, family = "Arial"),
      plot.margin      = margin(2, 4, 2, 2)
    )
}

# ── Panel A: Lac heatmap ─────────────────────────────────────────────────────
pA <- ggplot(lac_grid, aes(short, display, fill = frac)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = if_else(no_growth, "—",
                                 sprintf("%.0f", frac * 100))),
            size = 1.8,
            colour = if_else(lac_grid$no_growth, "grey55",
                             if_else(lac_grid$frac > 0.25, "white", "grey15"))) +
  facet_wrap(~ pool_label, nrow = 1) +
  scale_fill_gradient(low = "#f7f7f7", high = "#1b7837",
                      limits = c(0, NA),
                      labels = percent_format(accuracy = 1),
                      name = "Mutualism\n(% viable pairs)", na.value = "grey88") +
  labs(title = "Lactobacillus strain-resolved mutualism across the prebiotic gradient",
       subtitle = "— = strain non-viable on that medium",
       x = NULL, y = NULL) +
  theme_hm() +
  theme(legend.position = "right")

# ── Panel B: Per-strain trajectories ─────────────────────────────────────────
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

lac_gut_lines <- lac_grid |>
  filter(pool_label == "Gut (UHGG)", !is.na(frac))

pB <- ggplot(lac_gut_lines, aes(level_idx, frac, colour = group,
                                 group = display)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 1.0) +
  { if (has_ggrepel)
      ggrepel::geom_text_repel(
        data = lac_gut_lines |> filter(level_idx == max(level_idx)),
        aes(label = display), size = 1.8, hjust = 0, nudge_x = 0.3,
        direction = "y", segment.size = 0.2, max.overlaps = 15,
        family = "Arial")
    else
      geom_text(data = lac_gut_lines |> filter(level_idx == max(level_idx)),
                aes(label = display), size = 1.8, hjust = 0, nudge_x = 0.3)
  } +
  scale_colour_manual(values = GROUP_COLOURS, name = "Niche") +
  scale_x_continuous(breaks = LEVELS$level_idx, labels = LEVELS$short,
                     expand = expansion(mult = c(0.03, 0.25))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "Strain trajectories (gut partners)",
       x = NULL, y = "Mutualism %") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 7.5),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.15),
        legend.position = "left",
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        plot.tag = element_text(face = "bold", size = 9, family = "Arial"),
        plot.margin = margin(2, 4, 2, 2))

# ── Panel C: Akk heatmap ────────────────────────────────────────────────────
pC <- ggplot(akk_grid, aes(short, display, fill = frac)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = if_else(no_growth, "—",
                                 sprintf("%.0f", frac * 100))),
            size = 2.0,
            colour = if_else(akk_grid$no_growth, "grey55",
                             if_else(akk_grid$frac > 0.12, "white", "grey15"))) +
  scale_fill_gradient(low = "#f7f7f7", high = "#762a83",
                      limits = c(0, NA),
                      labels = percent_format(accuracy = 1),
                      name = "Mutualism\n(% viable pairs)", na.value = "grey88") +
  labs(title = "Akkermansia strain-resolved mutualism (gut partners)",
       x = NULL, y = NULL) +
  theme_hm() +
  theme(legend.position = "right")

# ── Assemble ─────────────────────────────────────────────────────────────────
fig6 <- pA / pC +
  plot_layout(heights = c(1.3, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 9, family = "Arial"))

ggsave(file.path(OUTDIR, "fig6_strain_heatmap.tiff"), fig6,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "fig6_strain_heatmap.pdf"), fig6,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "fig6_strain_heatmap.png"), fig6,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

# Summary table
summary_tbl <- lac_grid |>
  filter(!is.na(frac)) |>
  group_by(display, group, pool_label) |>
  summarise(mean_mut = round(mean(frac) * 100, 1),
            max_mut  = round(max(frac) * 100, 1),
            .groups = "drop") |>
  arrange(pool_label, desc(mean_mut))
write_tsv(summary_tbl, file.path(OUTDIR, "fig6_strain_summary.tsv"))

cat("Saved fig6_strain_heatmap (TIFF/PDF/PNG) to", OUTDIR, "\n")
