#!/usr/bin/env Rscript
# =============================================================================
# plot_fig4_mechanisms.R   (Main Figure 4)
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
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

THRESHOLD <- 1e-4
FIG_W <- 170 / 25.4
FIG_H <- 145 / 25.4
FIG_DPI <- 300

LEVELS <- tibble::tribble(
  ~level_idx, ~level_file,              ~short,
  0, "L0_base","L0", 1, "L1_inulin","L1", 2, "L2_fos","L2",
  3, "L3_gos","L3", 4, "L4_xos","L4", 5, "L5_pectin","L5",
  6, "L6_resistant_starch","L6", 7, "L7_bglucan","L7",
  8, "L8_hmo","L8", 9, "L9_mos","L9"
)

SYSTEMS <- tibble::tribble(
  ~sys_dir,        ~label,              ~colour,
  "akk_vs_uhgg",   "Akk × Gut",        "#762a83",
  "lac_vs_uhgg",   "Lac × Gut",         "#1b7837"
)
sys_colours <- setNames(SYSTEMS$colour, SYSTEMS$label)

# ── Load ─────────────────────────────────────────────────────────────────────
load_system <- function(sys_dir, label) {
  rows <- list(); frac_rows <- list()
  for (l in seq_len(nrow(LEVELS))) {
    p <- file.path(RESULTS, sys_dir, paste0(LEVELS$level_file[l], ".tsv"))
    if (!file.exists(p)) next
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    viable <- df |> filter(growth_a_alone > THRESHOLD, growth_b_alone > THRESHOLD)
    n_v <- nrow(viable); if (n_v == 0) next
    m <- viable |> summarise(
      ci_mean = mean(competition_intensity, na.rm = TRUE),
      ci_sd   = sd(competition_intensity, na.rm = TRUE),
      gs_mean = mean(gene_supported_fraction, na.rm = TRUE),
      gs_sd   = sd(gene_supported_fraction, na.rm = TRUE),
      ne_mean = mean(n_exchanged_metabolites, na.rm = TRUE),
      ne_sd   = sd(n_exchanged_metabolites, na.rm = TRUE),
      ga_mean = mean(growth_a_alone), gb_mean = mean(growth_b_alone),
      gca_mean = mean(growth_a_co), gcb_mean = mean(growth_b_co),
      ba_mean = mean(benefit_a), bb_mean = mean(benefit_b))
    rows[[length(rows)+1]] <- bind_cols(
      tibble(system = label, level_idx = LEVELS$level_idx[l], n = n_v), m)
    fc <- viable |> count(interaction_type) |>
      mutate(frac = n / sum(n)) |> select(interaction_type, frac) |>
      pivot_wider(names_from = interaction_type, values_from = frac, values_fill = 0)
    frac_rows[[length(frac_rows)+1]] <- bind_cols(
      tibble(system = label, level_idx = LEVELS$level_idx[l]), fc)
  }
  list(cont = bind_rows(rows), frac = bind_rows(frac_rows))
}

all_d <- map2(SYSTEMS$sys_dir, SYSTEMS$label, load_system)
cont <- bind_rows(map(all_d, "cont")) |>
  left_join(LEVELS, by = "level_idx") |> mutate(short = factor(short, levels = LEVELS$short))
frac <- bind_rows(map(all_d, "frac")) |>
  left_join(LEVELS, by = "level_idx") |> mutate(short = factor(short, levels = LEVELS$short))
if (!"mutualism" %in% names(frac)) frac$mutualism <- 0

# ── Theme ────────────────────────────────────────────────────────────────────
theme_gm <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title         = element_text(face = "bold", size = 7, hjust = 0,
                                        margin = margin(0, 0, 1, 0)),
      axis.title         = element_text(size = 6.5),
      axis.text          = element_text(size = 6, colour = "black"),
      axis.line          = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks         = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.15),
      legend.position    = "none",
      plot.tag           = element_text(face = "bold", size = 9, family = "Arial"),
      plot.margin        = margin(2, 4, 2, 2)
    )
}

x_sc <- scale_x_continuous(breaks = LEVELS$level_idx, labels = LEVELS$short,
                           expand = expansion(add = 0.3))
x_lab <- labs(x = "Prebiotic gradient (L0–L9)")
glucose_vline <- geom_vline(xintercept = 5.5, linetype = "dashed",
                            colour = "#CC3311", linewidth = 0.3)

# ── A: Competition intensity ────────────────────────────────────────────────
pA <- ggplot(cont, aes(level_idx, ci_mean, colour = system)) +
  geom_ribbon(aes(ymin = pmax(ci_mean - ci_sd, 0), ymax = ci_mean + ci_sd,
                  fill = system), alpha = 0.10, colour = NA) +
  glucose_vline +
  geom_line(linewidth = 0.5) + geom_point(size = 1.0) +
  scale_colour_manual(values = sys_colours) +
  scale_fill_manual(values = sys_colours) + x_sc + x_lab +
  labs(title = "Competition intensity", y = "C (resource overlap)") +
  theme_gm()

# ── B: Cross-fed metabolites ────────────────────────────────────────────────
pB <- ggplot(cont, aes(level_idx, ne_mean, colour = system)) +
  geom_ribbon(aes(ymin = pmax(ne_mean - ne_sd, 0), ymax = ne_mean + ne_sd,
                  fill = system), alpha = 0.10, colour = NA) +
  glucose_vline +
  geom_line(linewidth = 0.5) + geom_point(size = 1.0) +
  scale_colour_manual(values = sys_colours) +
  scale_fill_manual(values = sys_colours) + x_sc + x_lab +
  labs(title = "Cross-fed metabolites", y = "Metabolites / pair") +
  theme_gm()

# ── C: Gene-supported fraction ──────────────────────────────────────────────
pC <- ggplot(cont, aes(level_idx, gs_mean, colour = system)) +
  geom_ribbon(aes(ymin = pmax(gs_mean - gs_sd, 0),
                  ymax = pmin(gs_mean + gs_sd, 1),
                  fill = system), alpha = 0.10, colour = NA) +
  glucose_vline +
  geom_line(linewidth = 0.5) + geom_point(size = 1.0) +
  scale_colour_manual(values = sys_colours) +
  scale_fill_manual(values = sys_colours) + x_sc + x_lab +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Gene-supported fraction", y = "Fraction") +
  theme_gm()

# ── D: Growth ratio ─────────────────────────────────────────────────────────
ratio_df <- cont |>
  mutate(ra = gca_mean / pmax(ga_mean, 1e-6),
         rb = gcb_mean / pmax(gb_mean, 1e-6)) |>
  select(level_idx, short, system, ra, rb) |>
  pivot_longer(c(ra, rb), names_to = "sp", values_to = "ratio") |>
  mutate(sp_lab = if_else(sp == "ra", "Probiotic", "Commensal"))

pD <- ggplot(ratio_df, aes(level_idx, ratio, colour = system, linetype = sp_lab)) +
  geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.25) +
  glucose_vline +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
  scale_colour_manual(values = sys_colours) +
  scale_linetype_manual(values = c("Probiotic" = "solid", "Commensal" = "dashed"),
                        name = NULL) + x_sc + x_lab +
  labs(title = expression(paste(bold("Growth ratio  "), mu[co], "/", mu[alone])),
       y = expression(mu[co] / mu[alone])) +
  theme_gm() +
  theme(legend.position = c(0.78, 0.88),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ── E: Monoculture growth ───────────────────────────────────────────────────
gr_df <- cont |>
  select(level_idx, short, system, ga_mean, gb_mean) |>
  pivot_longer(c(ga_mean, gb_mean), names_to = "sp", values_to = "growth") |>
  mutate(sp_lab = if_else(sp == "ga_mean", "Probiotic", "Commensal"))

pE <- ggplot(gr_df, aes(level_idx, growth, colour = system, linetype = sp_lab)) +
  glucose_vline +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
  scale_colour_manual(values = sys_colours) +
  scale_linetype_manual(values = c("Probiotic" = "solid", "Commensal" = "dashed"),
                        name = NULL) + x_sc + x_lab +
  labs(title = "Monoculture growth rates",
       y = expression(paste(mu[alone], " (h"^-1, ")"))) +
  theme_gm() +
  theme(legend.position = c(0.22, 0.88),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ── F: Δ Mutualism vs Δ Competition ─────────────────────────────────────────
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

trans_df <- frac |>
  left_join(cont |> select(level_idx, system, ci_mean), by = c("level_idx", "system")) |>
  arrange(system, level_idx) |> group_by(system) |>
  mutate(dm = mutualism - lag(mutualism),
         dc = ci_mean - lag(ci_mean),
         tr = paste0(lag(short), "→", short)) |>
  filter(!is.na(dm)) |> ungroup()
crash <- trans_df |> filter(level_idx == 6)

xr <- range(trans_df$dc, na.rm = TRUE); yr <- range(trans_df$dm, na.rm = TRUE)
pF <- ggplot(trans_df, aes(dc, dm, colour = system)) +
  # shade the antagonistic quadrant (more competition, less mutualism) where the crash lands
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#CC3311", alpha = 0.06) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.2) +
  # directional arrow from the cluster toward the crash corner
  annotate("segment", x = xr[2]*0.18, y = yr[1]*0.10,
           xend = xr[2]*0.80, yend = yr[1]*0.82,
           colour = "#CC3311", linewidth = 0.3, alpha = 0.5,
           arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  annotate("text", x = xr[2]*0.50, y = yr[1]*0.20,
           label = "more competition,\nless cooperation",
           size = 1.8, colour = "#CC3311", fontface = "italic",
           lineheight = 0.85, hjust = 0.5) +
  geom_point(size = 1.4, alpha = 0.75) +
  geom_point(data = crash, size = 3.0, shape = 1, stroke = 0.8) +
  { if (has_ggrepel)
      ggrepel::geom_text_repel(data = crash, aes(label = "L5→L6"),
                                size = 2, fontface = "bold",
                                colour = "#CC3311", family = "Arial",
                                nudge_x = 0.0015, nudge_y = -0.004,
                                min.segment.length = 0)
    else
      geom_text(data = crash, aes(label = "L5→L6"),
                size = 2, fontface = "bold", vjust = -1, colour = "#CC3311")
  } +
  scale_colour_manual(values = sys_colours, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(title = expression(paste(Delta, bold(" Mutualism vs "), Delta, bold(" Competition"))),
       subtitle = "Each point = one gradient step (L→L+1); rings = L5→L6 (free glucose enters)",
       x = expression(Delta ~ C ~ "(competition intensity)"),
       y = expression(Delta ~ "Mutualism fraction")) +
  theme_gm() +
  theme(plot.subtitle = element_text(size = 5, colour = "grey35"),
        legend.position = c(0.80, 0.82),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ── Assemble ─────────────────────────────────────────────────────────────────
# Panel order follows the order panels are first cited in the text:
#   A Competition intensity · B Δ Mutualism vs Δ Competition (glucose crash) ·
#   C Monoculture growth rates · D Growth ratio μco/μalone ·
#   E Cross-fed metabolites · F Gene-supported fraction.
#   tag_levels auto-labels in this composition order.
fig4 <- (pA + pF + pE) / (pD + pB + pC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 9, family = "Arial"))

ggsave(file.path(OUTDIR, "fig4_mechanisms.tiff"), fig4,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "fig4_mechanisms.pdf"), fig4,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "fig4_mechanisms.png"), fig4,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved fig4_mechanisms (TIFF/PDF/PNG) to", OUTDIR, "\n")
