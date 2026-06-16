#!/usr/bin/env Rscript
# =============================================================================
# plot_fig3_gradient_overview.R   (Main Figure 3)
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
FIG_H <- 160 / 25.4
FIG_DPI <- 300

LEVELS <- tibble::tribble(
  ~level_idx, ~level_file,              ~short,     ~display,       ~category,
  0, "L0_base",             "L0",  "L0\nBase",       "baseline",
  1, "L1_inulin",           "L1",  "L1\nInulin",     "universal",
  2, "L2_fos",              "L2",  "L2\nFOS",        "universal",
  3, "L3_gos",              "L3",  "L3\nGOS",        "specialty",
  4, "L4_xos",              "L4",  "L4\nXOS",        "specialty",
  5, "L5_pectin",           "L5",  "L5\nPectin",     "specialty",
  6, "L6_resistant_starch", "L6",  "L6\nRS",         "glucose",
  7, "L7_bglucan",          "L7",  "L7\nβ-glucan","saturated",
  8, "L8_hmo",              "L8",  "L8\nHMO",        "saturated",
  9, "L9_mos",              "L9",  "L9\nMOS",        "saturated"
)

SYSTEMS <- tibble::tribble(
  ~sys_dir,        ~label,                            ~colour,
  "akk_vs_uhgg",   "Akkermansia × Gut (UHGG)",       "#762a83",
  "lac_vs_uhgg",   "Lactobacillus × Gut (UHGG)",     "#1b7837"
)
sys_colours <- setNames(SYSTEMS$colour, SYSTEMS$label)

PHASE_COLOURS <- c(baseline = "#999999", universal = "#D55E00",
                   specialty = "#009E73", glucose = "#CC3311",
                   saturated = "#56B4E9")

# ── Load ─────────────────────────────────────────────────────────────────────
load_system <- function(sys_dir, label) {
  rows <- list()
  for (l in seq_len(nrow(LEVELS))) {
    p <- file.path(RESULTS, sys_dir, paste0(LEVELS$level_file[l], ".tsv"))
    if (!file.exists(p)) { warning("Missing: ", p); next }
    df <- suppressMessages(read_tsv(p, show_col_types = FALSE, progress = FALSE))
    viable <- df |> filter(growth_a_alone > THRESHOLD, growth_b_alone > THRESHOLD)
    n_viable <- nrow(viable)
    if (n_viable == 0) next
    counts <- viable |> count(interaction_type) |>
      mutate(frac = n / sum(n)) |> select(interaction_type, frac) |>
      pivot_wider(names_from = interaction_type, values_from = frac, values_fill = 0)
    metrics <- viable |> summarise(
      benefit_a_mean = mean(benefit_a, na.rm = TRUE),
      benefit_b_mean = mean(benefit_b, na.rm = TRUE))
    rows[[length(rows)+1]] <- bind_cols(
      tibble(system = label, level_idx = LEVELS$level_idx[l], n_viable = n_viable),
      counts, metrics)
  }
  bind_rows(rows)
}

all_summ <- bind_rows(map2(SYSTEMS$sys_dir, SYSTEMS$label, load_system)) |>
  left_join(LEVELS, by = "level_idx") |>
  mutate(display = factor(display, levels = LEVELS$display))
if (!"mutualism"  %in% names(all_summ)) all_summ$mutualism  <- 0
if (!"competition" %in% names(all_summ)) all_summ$competition <- 0
all_summ <- all_summ |> mutate(net_benefit = (benefit_a_mean + benefit_b_mean) / 2)
cat("Loaded", nrow(all_summ), "system × level rows\n")

# ── Theme ────────────────────────────────────────────────────────────────────
theme_gm <- function() {
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      plot.title         = element_text(face = "bold", size = 8, hjust = 0,
                                        margin = margin(0, 0, 2, 0)),
      axis.title         = element_text(size = 7),
      axis.text          = element_text(size = 6, colour = "black"),
      axis.text.x        = element_text(lineheight = 0.85),
      axis.line          = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks         = element_line(linewidth = 0.3, colour = "black"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 6),
      legend.key.size    = unit(0.3, "cm"),
      legend.key.height  = unit(0.3, "cm"),
      legend.margin      = margin(0, 0, 0, 0),
      plot.tag           = element_text(face = "bold", size = 10),
      plot.margin        = margin(3, 5, 3, 3)
    )
}

band_df <- LEVELS |> mutate(xmin = level_idx - 0.5, xmax = level_idx + 0.5)
phase_bg <- function() {
  list(
    geom_rect(data = band_df,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category),
              alpha = 0.10, inherit.aes = FALSE),
    scale_fill_manual(values = PHASE_COLOURS, guide = "none")
  )
}

# Glucose vertical line at L5.5 — shared across all panels
glucose_vline <- geom_vline(xintercept = 5.5, linetype = "dashed",
                            colour = "#CC3311", linewidth = 0.4)

x_scale <- scale_x_continuous(breaks = LEVELS$level_idx, labels = LEVELS$short,
                              expand = expansion(add = 0.5))
x_lab <- labs(x = "Prebiotic gradient (L0–L9)")

# ── Panel A: Viability ───────────────────────────────────────────────────────
pA <- ggplot(all_summ, aes(level_idx, n_viable, colour = system)) +
  phase_bg() +
  glucose_vline +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  # Glucose annotation
  annotate("text", x = 5.8, y = max(all_summ$n_viable) * 0.55,
           label = "Glucose\nadded", size = 2, colour = "#CC3311",
           fontface = "bold", lineheight = 0.85, hjust = 0, family = "Arial") +
  scale_colour_manual(values = sys_colours) +
  x_scale + x_lab +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Viable pair count across gradient", y = "Viable pairs") +
  theme_gm() +
  theme(legend.position = c(0.60, 0.30),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ── Panel B: Mutualism fraction ──────────────────────────────────────────────
# Find peak per system
mut_peaks <- all_summ |> group_by(system) |>
  slice_max(mutualism, n = 1) |> ungroup()

pB <- ggplot(all_summ, aes(level_idx, mutualism, colour = system)) +
  phase_bg() +
  glucose_vline +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  geom_point(data = mut_peaks, size = 2.5, shape = 1, stroke = 0.8) +
  # Peak annotation at top
  annotate("label", x = 5, y = max(all_summ$mutualism, na.rm = TRUE) * 1.12,
           label = "Peak at L5 (Pectin)", size = 2, colour = "#009E73",
           fontface = "bold", family = "Arial",
           fill = alpha("white", 0.9), label.size = 0.2,
           label.padding = unit(0.12, "lines")) +
  scale_colour_manual(values = sys_colours) +
  x_scale + x_lab +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = "Mutualism peaks at specialty sugars",
       y = "Mutualism (% of viable pairs)") +
  theme_gm() +
  theme(legend.position = "none")

# ── Panel C: Competition fraction ────────────────────────────────────────────
pC <- ggplot(all_summ, aes(level_idx, competition, colour = system)) +
  phase_bg() +
  glucose_vline +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = sys_colours) +
  x_scale + x_lab +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.03, 0.03))) +
  labs(title = "Competition mirrors mutualism decline",
       y = "Competition (% of viable pairs)") +
  theme_gm() +
  theme(legend.position = "none")

# ── Panel D: Mean net benefit ────────────────────────────────────────────────
pD <- ggplot(all_summ, aes(level_idx, net_benefit, colour = system)) +
  phase_bg() +
  glucose_vline +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, linetype = "solid") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = sys_colours) +
  x_scale + x_lab +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Net benefit least negative at specialty sugars",
       y = expression(paste("Mean  ", bar(beta), " = (", beta[A], "+", beta[B], ")/2"))) +
  theme_gm() +
  theme(legend.position = "none")

# ── Assemble ─────────────────────────────────────────────────────────────────
fig3 <- (pA + pB) / (pC + pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10, family = "Arial"))

ggsave(file.path(OUTDIR, "fig3_gradient_overview.tiff"), fig3,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       compression = "lzw", device = "tiff")
ggsave(file.path(OUTDIR, "fig3_gradient_overview.pdf"), fig3,
       width = FIG_W, height = FIG_H, device = cairo_pdf)
ggsave(file.path(OUTDIR, "fig3_gradient_overview.png"), fig3,
       width = FIG_W, height = FIG_H, dpi = FIG_DPI)

cat("Saved fig3_gradient_overview (TIFF/PDF/PNG) to", OUTDIR, "\n")
