#!/usr/bin/env Rscript
# =============================================================================
# plot_fig4_mechanisms.R   (Main Figure 4)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
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


args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { i <- match(f, args); if (!is.na(i) && i < length(args)) args[i+1] else d }
RESULTS <- get_arg("--results", "results")
OUTDIR  <- get_arg("--outdir",  file.path(RESULTS, "figures_paper"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

THRESHOLD <- 1e-4
FIG_W <- 170 / 25.4
FIG_H <- 225 / 25.4
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
  "akk_vs_uhgg",   "Akkermansia × Gut",   "#762a83",
  "lac_vs_uhgg",   "Lactobacillus × Gut", "#1b7837"
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
  theme_classic(base_size = 7, base_family = BMC_FONT) +
    theme(
      plot.title         = md(face = "bold", size = 7, hjust = 0,
                                        margin = margin(0, 0, 1, 0)),
      axis.title         = md(size = 6.5),
      axis.title.x       = md(size = 6.5), axis.title.y = md(size = 6.5, angle = 90),
      axis.text          = element_text(size = 6, colour = "black"),
      axis.line          = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks         = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.15),
      legend.text        = md(),
      legend.position    = "none",
      plot.tag           = element_text(face = "bold", size = 9, family = BMC_FONT),
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
  scale_colour_manual(values = sys_colours, labels = sp_md, name = NULL) +
  scale_linetype_manual(values = c("Probiotic" = "solid", "Commensal" = "dashed"),
                        name = NULL) + x_sc + x_lab +
  labs(title = "**Growth ratio** μ<sub>co</sub>/μ<sub>alone</sub>",
       y = "μ<sub>co</sub>/μ<sub>alone</sub>") +
  theme_gm() +
  theme(legend.position = c(0.80, 0.64), legend.spacing.y = unit(0.02, "cm"),
        legend.text = md(size = 5.5),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ── E: Monoculture growth ───────────────────────────────────────────────────
gr_df <- cont |>
  select(level_idx, short, system, ga_mean, gb_mean)
# The partner curve is the same in both systems: every probiotic strain is viable at
# every level, so both systems average the same set of viable UHGG partners (each
# weighted equally). Draw it once, in grey, instead of hiding one line under the other.
gb_wide <- tidyr::pivot_wider(gr_df |> select(level_idx, system, gb_mean),
                              names_from = system, values_from = gb_mean)
stopifnot(max(abs(gb_wide[[2]] - gb_wide[[3]])) < 1e-9)
PARTNER <- "Gut partners (both systems)"
gr_long <- bind_rows(
  gr_df |> transmute(level_idx, growth = ga_mean,
                     series = paste0(sub(" × Gut", "", system), " (probiotic)")),
  gr_df |> filter(system == SYSTEMS$label[1]) |>
    transmute(level_idx, growth = gb_mean, series = PARTNER)) |>
  mutate(series = factor(series, levels = c("Akkermansia (probiotic)",
                                            "Lactobacillus (probiotic)", PARTNER)))
ser_col <- setNames(c(SYSTEMS$colour, "grey45"), levels(gr_long$series))
ser_lt  <- setNames(c("solid", "solid", "dashed"), levels(gr_long$series))

pE <- ggplot(gr_long, aes(level_idx, growth, colour = series, linetype = series)) +
  glucose_vline +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
  scale_colour_manual(values = ser_col, labels = sp_md, name = NULL) +
  scale_linetype_manual(values = ser_lt, labels = sp_md, name = NULL) + x_sc + x_lab +
  labs(title = "Monoculture growth rates",
       y = "μ<sub>alone</sub> (h<sup>−1</sup>)") +
  theme_gm() +
  theme(legend.position = c(1, 0.02), legend.justification = c(1, 0),
        legend.text = md(size = 5.5),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_blank(), legend.key = element_blank())

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
           size = 1.8, colour = "#CC3311",
           lineheight = 0.85, hjust = 0.5) +
  geom_point(size = 1.4, alpha = 0.75) +
  geom_point(data = crash, size = 3.0, shape = 1, stroke = 0.8) +
  { if (has_ggrepel)
      ggrepel::geom_text_repel(data = crash, aes(label = "L5→L6"),
                                size = 2, fontface = "bold",
                                colour = "#CC3311", family = BMC_FONT,
                                nudge_x = 0.0015, nudge_y = -0.004,
                                min.segment.length = 0, seed = 1)
    else
      geom_text(data = crash, aes(label = "L5→L6"),
                size = 2, fontface = "bold", vjust = -1, colour = "#CC3311")
  } +
  scale_colour_manual(values = sys_colours, name = NULL, labels = sp_md) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(title = "Δ **Mutualism vs** Δ **Competition**",
       subtitle = "Each point = one gradient step (L→L+1); rings = L5→L6 (free glucose enters)",
       x = "Δ C (competition intensity)",
       y = "Δ Mutualism fraction") +
  theme_gm() +
  theme(plot.subtitle = element_text(size = 5, colour = "grey35"),
        legend.position = c(0.80, 0.82),
        legend.text = md(size = 5.5),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ── Assemble ─────────────────────────────────────────────────────────────────
# ── Panel G: controlled carbon-quality contrast (elevated from former Fig. S7) ──
# Same L1–L4 background; only the last carbon changes (pectin L5 -> glucose L5-glc),
# holding richness constant, so this isolates carbon QUALITY from the L5->L6 step
# (which also adds maltose + maltodextrin). Both systems drop.
gi_spec <- tibble::tribble(~system,~pre, "Akkermansia × Gut","akk", "Lactobacillus × Gut","lac")
gi_lvl  <- tibble::tribble(~lvl,~file,~cls,
  "L5\n(pectin)","L5_pectin","complex carbon",
  "L5-glc\n(glucose)","L5glc","free glucose",
  "L6\n(+free sugars)","L6_resistant_starch","more free sugars")
gi_mut <- function(pre, file) {
  d <- suppressMessages(read_tsv(file.path(RESULTS, "fig4",
            paste0(pre, "_", file, ".tsv")), show_col_types = FALSE, progress = FALSE)) |>
    filter(growth_a_alone > THRESHOLD, growth_b_alone > THRESHOLD)
  100 * mean(d$interaction_type == "mutualism")
}
gi <- tidyr::crossing(gi_spec, gi_lvl) |> rowwise() |>
  mutate(mut = gi_mut(pre, file)) |> ungroup() |>
  mutate(system = factor(system, levels = c("Akkermansia × Gut","Lactobacillus × Gut")),
         lvl = factor(lvl, levels = unique(gi_lvl$lvl)),
         cls = factor(cls, levels = c("complex carbon","free glucose","more free sugars")))
gi_cols <- c("complex carbon"="#41ab5d","free glucose"="#e6862e","more free sugars"="#d6604d")
pG <- ggplot(gi, aes(lvl, mut, fill = cls)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", mut)), vjust = -0.3, size = 1.9) +
  facet_wrap(~system, scales = "free_y", nrow = 1, labeller = sp_labeller) +
  scale_fill_manual(values = gi_cols, name = "Last carbon") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Controlled contrast (richness fixed): pectin → glucose lowers mutualism in both systems",
       x = NULL, y = "Mutualism (%)") +
  theme_gm() +
  theme(legend.position = "right", legend.key.size = unit(0.22, "cm"),
        legend.text = md(size = 5.5), legend.title = element_text(size = 6),
        strip.text = md(face = "bold", size = 6.5),
        axis.text.x = element_text(size = 5, lineheight = 0.8))

# Panel order follows the order panels are first cited in the text:
#   A Competition intensity · B Δ Mutualism vs Δ Competition (glucose crash) ·
#   C Monoculture growth rates · D Growth ratio μco/μalone ·
#   E Cross-fed metabolites · F Gene-supported fraction · G controlled pectin→glucose contrast.
#   Layout: colour key on top, then two panels per row (A B / C D / E F / G).

# Figure-level colour key: panels A, E and F carry no legend of their own.
# Drawn as a small plot with fixed positions (swatch = ribbon + line + point),
# so no legend text has to be measured.
key_df <- tibble::tibble(x0 = c(0.15, 2.55), lab_x = c(0.75, 3.15),
                         colour = unname(sys_colours),
                         lab = c("*Akkermansia* \u00d7 Gut", "*Lactobacillus* \u00d7 Gut"))
key_row <- ggplot(key_df) +
  geom_rect(aes(xmin = x0, xmax = x0 + 0.5, ymin = 0.3, ymax = 0.7, fill = colour), alpha = 0.25) +
  geom_segment(aes(x = x0, xend = x0 + 0.5, y = 0.5, yend = 0.5, colour = colour), linewidth = 0.5) +
  geom_point(aes(x = x0 + 0.25, y = 0.5, colour = colour), size = 1.0) +
  ggtext::geom_richtext(aes(x = lab_x, y = 0.5, label = lab), hjust = 0, size = 6.5 / .pt,
                        family = BMC_FONT, fill = NA, label.colour = NA,
                        label.padding = unit(0, "pt")) +
  annotate("text", x = 10, y = 0.5, hjust = 1, size = 5.5 / .pt, colour = "grey35",
           family = BMC_FONT,
           label = "line = mean across viable pairs; shading (A, E, F) = \u00b1 SD") +
  scale_colour_identity() + scale_fill_identity() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void(base_family = BMC_FONT)

# Tags are set by hand so that the colour key is not lettered.
tagged <- function(p, t) p + labs(tag = t) +
  theme(plot.tag = element_text(face = "bold", size = 9, family = BMC_FONT))
fig4 <- key_row / (tagged(pA, "A") + tagged(pF, "B")) / (tagged(pE, "C") + tagged(pD, "D")) /
  (tagged(pB, "E") + tagged(pC, "F")) / tagged(pG, "G") +
  plot_layout(heights = c(0.07, 1, 1, 1, 0.85))

ggsave(file.path(OUTDIR, "fig4_mechanisms.tiff"), fig4,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "fig4_mechanisms.pdf"), fig4,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "fig4_mechanisms.png"), fig4,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved fig4_mechanisms (TIFF/PDF/PNG) to", OUTDIR, "\n")
