#!/usr/bin/env Rscript
# Combined benchmark figure: speed scaling + memory + time breakdown + correctness.
#
# Final 2×2 layout (panel tags follow text citation order):
#   A  Correctness scatter (numerical agreement)        [built as pE]
#   B  Wall time vs threads (log10 y) + speedup labels   [built as pA]
#   C  Speedup vs threads — fast-mic only                [built as pB]
#   D  Peak memory vs threads (+ COBRApy 1-thread ref)   [built as pC]
# (pD, the load-vs-FBA time breakdown, is built but not placed in the figure.)
#
# Usage:
#   Rscript plot_benchmark.R \
#     results/fig1/benchmark/thread_scaling/thread_scaling_results.tsv \
#     results/fig1/benchmark/correctness/scatter.tsv \
#     results/fig1/benchmark/correctness/stats.tsv \
#     [output_prefix]
#
#   output_prefix defaults to "results/figures_paper/fig1_benchmark"

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggrepel)
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


# ── CLI args ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
scaling_tsv  <- if (length(args) >= 1) args[1] else
  "results/fig1/benchmark/thread_scaling/thread_scaling_results.tsv"
scatter_tsv  <- if (length(args) >= 2) args[2] else
  "results/fig1/benchmark/correctness/scatter.tsv"
stats_tsv    <- if (length(args) >= 3) args[3] else
  "results/fig1/benchmark/correctness/stats.tsv"
out_prefix   <- if (length(args) >= 4) args[4] else
  "results/figures_paper/fig1_benchmark"

# ── Shared theme ──────────────────────────────────────────────────────────────
FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 225 / 25.4      # BMC maximum height

base_theme <- theme_bw(base_size = 7, base_family = BMC_FONT) +
  theme(
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom",
    legend.key.size   = unit(0.28, "cm"),
    legend.text       = md(size = 5.8),
    legend.title      = element_text(size = 6.3),
    axis.text         = element_text(size = 5.5),
    plot.title        = element_text(size = 7.5, face = "bold"),
    plot.subtitle     = element_text(size = 5.8, colour = "grey35"),
    axis.title        = md(size = 6.5),
    axis.title.x      = md(size = 6.5), axis.title.y = md(size = 6.5, angle = 90),
    strip.text        = md(size = 6.5),
    plot.tag          = element_text(size = 9, face = "bold")
  )

BLUE   <- "#2166AC"
RED    <- "#D6604D"
GREY   <- "#999999"
GREEN  <- "#4DAC26"
ORANGE <- "#E08214"

# ── Helper: read stat value ───────────────────────────────────────────────────
get_stat <- function(stats, key) {
  v <- stats$value[stats$metric == key]
  if (!length(v)) return(NA_real_)
  suppressWarnings(as.numeric(v))
}

# ============================================================
# Load data
# ============================================================

scaling <- tryCatch(read.delim(scaling_tsv, stringsAsFactors = FALSE),
                    error = function(e) NULL)
scatter  <- tryCatch(read.delim(scatter_tsv, stringsAsFactors = FALSE,
                                na.strings = "NA"),
                     error = function(e) NULL)
stats    <- tryCatch(read.delim(stats_tsv,   stringsAsFactors = FALSE),
                     error = function(e) NULL)

if (is.null(scaling) || nrow(scaling) == 0)
  stop("No thread-scaling data found in: ", scaling_tsv)
if (is.null(scatter) || nrow(scatter) == 0)
  stop("No correctness scatter data found in: ", scatter_tsv)

# Detect available columns — older TSVs may lack avg_load_s / avg_fba_s
has_timing_cols <- all(c("avg_load_s", "avg_fba_s") %in% names(scaling))

# ── Numeric n_models ordering — applied EVERYWHERE so legends and facets
#    don't fall into character sort ("100" < "1000" < "500"). ──────────────
n_sizes_int <- sort(unique(as.integer(scaling$n_models)))
n_labels    <- as.character(n_sizes_int)
make_nfac   <- function(x) factor(as.character(x), levels = n_labels)

# How many repeats per (tool, n, threads)?  Used in the caption.
n_repeats <- scaling |>
  group_by(tool, n_models, threads) |>
  summarise(k = n(), .groups = "drop") |>
  pull(k) |>
  median()

# ============================================================
# PANELS A, B, C  — thread scaling
# ============================================================

# Summarise per (tool, n_models, threads): min/median/max over repeats
all_summ <- scaling |>
  group_by(tool, n_models, threads) |>
  summarise(
    wall_min = min(wall_s),
    wall_med = median(wall_s),
    wall_max = max(wall_s),
    mem_min  = min(peak_mem_mb),
    mem_med  = median(peak_mem_mb),
    mem_max  = max(peak_mem_mb),
    .groups  = "drop"
  ) |>
  mutate(
    n_label  = make_nfac(n_models),
    grp      = paste(tool, n_label, sep = " · ")
  )

fm_summ <- all_summ |> filter(tool == "fast-mic")
cb_summ <- all_summ |> filter(tool == "COBRApy")

# Tool palette
tool_col <- c("fast-mic" = BLUE, "COBRApy" = RED)
size_lty <- setNames(
  rep(c("solid", "dashed", "dotdash", "longdash"), length.out = length(n_labels)),
  n_labels
)

x_breaks <- sort(unique(all_summ$threads))

# ── Panel A: Wall time (log10) — both tools + per-(n) speedup annotation ──
# Compute speedup of fast-mic over COBRApy at 1 thread, per n_models.
speedup_1t <- all_summ |>
  filter(threads == 1) |>
  select(tool, n_models, wall_med) |>
  pivot_wider(names_from = tool, values_from = wall_med) |>
  rename(fastmic = `fast-mic`, cobra = COBRApy) |>
  mutate(
    n_label  = make_nfac(n_models),
    speedup  = cobra / fastmic,
    label    = sprintf("%.0f× faster (n=%d)", speedup, n_models),
    # Midpoint between the two lines on log scale → geometric mean
    y_mid    = sqrt(cobra * fastmic)
  ) |>
  filter(!is.na(speedup))

pA <- ggplot(all_summ,
             aes(threads, wall_med,
                 colour = tool, linetype = n_label, group = grp)) +
  geom_ribbon(aes(ymin = wall_min, ymax = wall_max, fill = tool),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = wall_min, ymax = wall_max), width = 0.16,
                linewidth = 0.6, alpha = 0.95) +
  geom_point(size = 1.3) +
  {
    if (nrow(speedup_1t) > 0)
      geom_text(data = speedup_1t,
                aes(x = 1.05, y = y_mid, label = label),
                inherit.aes = FALSE,
                hjust = 0, vjust = -0.5, size = 1.7,
                colour = "grey20")
  } +
  scale_x_continuous(breaks = x_breaks, trans = "log2",
                     labels = function(x) as.integer(x)) +
  scale_y_log10(labels = label_number(suffix = " s", accuracy = 0.01)) +
  scale_colour_manual(values = tool_col, name = "Tool") +
  scale_fill_manual(values = tool_col, name = "Tool") +
  scale_linetype_manual(values = size_lty, name = "Models") +
  labs(x = "Threads", y = "Wall time (s, log scale)",
       title = "Runtime",
       subtitle = "fast-mic vs COBRApy at matched core counts") +
  base_theme

# ── Panel B: Speedup — fast-mic only ──
# COBRApy parallelism here comes from an external multiprocessing.Pool
# wrapper applied by the benchmark harness; COBRApy itself is single-threaded
# per FBA call.  Plotting its "speedup" would mislead readers, so we drop it.
fm_base <- fm_summ |> filter(threads == 1) |>
  select(n_models, base_wall = wall_med)

speedup_df <- fm_summ |>
  left_join(fm_base, by = "n_models") |>
  mutate(
    speedup_med = base_wall / wall_med,
    speedup_max = base_wall / wall_min,   # min wall → max speedup
    speedup_min = base_wall / wall_max,
    n_label     = make_nfac(n_models)
  )

ideal_df <- tibble(threads = x_breaks, speedup = as.numeric(x_breaks))

pB <- ggplot(speedup_df,
             aes(threads, speedup_med,
                 linetype = n_label, group = n_label)) +
  geom_ribbon(aes(ymin = speedup_min, ymax = speedup_max),
              fill = BLUE, alpha = 0.15, colour = NA) +
  geom_line(data = ideal_df, aes(threads, speedup),
            colour = GREY, linetype = "dashed", linewidth = 0.7,
            inherit.aes = FALSE) +
  annotate("text", x = max(x_breaks) * 0.9,
           y = max(x_breaks) * 0.92,
           label = "Ideal", colour = GREY, hjust = 1, size = 3) +
  geom_line(linewidth = 0.8, colour = BLUE) +
  geom_errorbar(aes(ymin = speedup_min, ymax = speedup_max), width = 0.16,
                linewidth = 0.6, colour = BLUE, alpha = 0.95) +
  geom_point(size = 1.3, colour = BLUE) +
  scale_x_continuous(breaks = x_breaks, trans = "log2",
                     labels = function(x) as.integer(x)) +
  scale_linetype_manual(values = size_lty, name = "Models") +
  labs(x = "Threads",
       y = "Speedup (×, relative to 1 thread)",
       title = "fast-mic parallel speedup",
       subtitle = "COBRApy omitted (single-threaded per FBA call)") +
  base_theme

# ── Panel C: Peak memory — fast-mic line + COBRApy 1-thread reference ──
# COBRApy 2/4/8/12-thread memory values are measurement artifacts of our
# multiprocessing wrapper (fork copy-on-write, lazy imports), not a
# property of COBRApy itself.  We keep only the 1-thread point as the
# honest single-process baseline.
mem_fm <- fm_summ   |> filter(mem_med > 0)
mem_cb <- cb_summ   |> filter(mem_med > 0, threads == 1) |>
  mutate(n_label = make_nfac(n_models))
# Label each fast-mic curve with its corpus size at the right-hand end,
# where the three n levels are well separated on the y-axis.
mem_fm_end <- mem_fm |>
  group_by(n_label) |>
  slice_max(threads, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(lbl = paste0("n=", as.character(n_label)))

pC <- ggplot(mem_fm,
             aes(threads, mem_med,
                 linetype = n_label, group = n_label)) +
  geom_ribbon(aes(ymin = mem_min, ymax = mem_max),
              fill = BLUE, alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8, colour = BLUE) +
  geom_errorbar(aes(ymin = mem_min, ymax = mem_max), width = 0.16,
                linewidth = 0.6, colour = BLUE, alpha = 0.95) +
  geom_point(size = 1.3, colour = BLUE) +
  geom_text(data = mem_fm_end,
            aes(x = threads, y = mem_med, label = lbl),
            inherit.aes = FALSE, colour = BLUE, size = 2.3,
            hjust = 0, nudge_x = 0.05, vjust = 0.4) +
  # COBRApy single-process reference points (at threads=1)
  {
    if (nrow(mem_cb) > 0)
      geom_point(data = mem_cb,
                 aes(x = threads, y = mem_med, shape = n_label),
                 colour = RED, size = 1.7, stroke = 0.6,
                 inherit.aes = FALSE)
  } +
  {
    if (nrow(mem_cb) > 0)
      geom_text_repel(data = mem_cb,
                aes(x = threads, y = mem_med,
                    label = sprintf("COBRApy 1t (n=%d)", n_models)),
                inherit.aes = FALSE, colour = RED, size = 2.3,
                direction = "y", hjust = 0, nudge_x = 0.3,
                box.padding = 0.25, min.segment.length = 0,
                segment.size = 0.3, segment.colour = RED, seed = 1)
  } +
  scale_x_continuous(breaks = x_breaks, trans = "log2",
                     labels = function(x) as.integer(x),
                     expand = expansion(mult = c(0.05, 0.16))) +
  scale_y_continuous(labels = label_number(suffix = " MB", accuracy = 1)) +
  scale_linetype_manual(values = size_lty, name = "Models") +
  scale_shape_manual(values = c(15, 17, 18, 19)[seq_along(n_labels)],
                     name = "Models") +
  labs(x = "Threads", y = "Peak memory (MB)",
       title = "Memory usage",
       subtitle = "fast-mic curve; COBRApy shown only at 1 thread (true baseline)") +
  base_theme

# ============================================================
# PANEL D — Stacked bar: load vs FBA time per thread, both tools
# ============================================================
#
# facet_grid(tool ~ n_models, scales = "free_y") gives every TOOL its own
# y-range (because COBRApy is ~50–100× slower than fast-mic) while all
# n_models columns within a tool share the same y — that way the reader
# can compare different model-set sizes without optical-illusion rescaling.

if (has_timing_cols) {
  bar_df <- scaling |>
    group_by(tool, n_models, threads) |>
    summarise(
      avg_load_s = mean(avg_load_s, na.rm = TRUE),
      avg_fba_s  = mean(avg_fba_s,  na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(cols = c(avg_load_s, avg_fba_s),
                 names_to = "component", values_to = "time_s") |>
    mutate(
      component = recode(component,
                         "avg_load_s" = "Model loading",
                         "avg_fba_s"  = "FBA solve"),
      n_label   = make_nfac(n_models),
      threads_f = factor(threads, levels = sort(unique(threads))),
      tool      = factor(tool, levels = c("fast-mic", "COBRApy"))
    ) |>
    filter(time_s > 0)

  pD <- ggplot(bar_df,
               aes(threads_f, time_s, fill = component)) +
    geom_col(position = "stack", width = 0.7) +
    facet_grid(tool ~ n_label,
               labeller = labeller(n_label = function(x) paste0("n=", x)),
               scales = "free_y") +
    scale_fill_manual(values = c("Model loading" = ORANGE, "FBA solve" = BLUE),
                      name = NULL) +
    scale_y_continuous(labels = label_number(suffix = " s", accuracy = 0.01)) +
    labs(x = "Threads", y = "Time per model (s)",
         title = "Per-model time breakdown",
         subtitle = "Same y-scale within each tool row; rows differ ~100×") +
    base_theme +
    theme(legend.position = "bottom",
          strip.text = md(size = 6))
} else {
  pD <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Timing breakdown\nnot available\n(re-run benchmark)", size = 4) +
    theme_void(base_family = BMC_FONT) +
    labs(title = "Time breakdown")
}

# ============================================================
# PANEL E — Correctness scatter
# ============================================================

r2   <- if (!is.null(stats)) get_stat(stats, "r_squared")  else NA
mae  <- if (!is.null(stats)) get_stat(stats, "mae")        else NA
n    <- if (!is.null(stats)) get_stat(stats, "n_models")   else nrow(scatter)
pct1 <- if (!is.null(stats)) get_stat(stats, "pct_agree_1pct") else NA

lim <- max(scatter$fastmic_growth, scatter$cobra_growth, na.rm = TRUE) * 1.05
if (!is.finite(lim) || lim == 0) lim <- 1

ann <- if (!is.na(r2) && !is.na(mae)) {
  sprintf("R2 = %.4f\nMAE = %.2e h-1\n%.1f%% within 1%%\nn = %d",
          r2, mae, pct1, as.integer(n))
} else {
  "No stats available"
}

# alpha=0.3 keeps overplotted dense regions readable at large n.
pE <- ggplot(scatter, aes(cobra_growth, fastmic_growth)) +
  geom_abline(slope = 1, intercept = 0,
              colour = GREY, linetype = "dashed", linewidth = 0.6) +
  geom_point(alpha = 0.30, size = 1.1, colour = BLUE) +
  annotate("text",
           x = lim * 0.04, y = lim * 0.97,
           label = ann, hjust = 0, vjust = 1,
           size = 2, family = BMC_FONT) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = "COBRApy growth rate (h<sup>−1</sup>)",
       y = "fast-mic growth rate (h<sup>−1</sup>)",
       title = "Numerical agreement") +
  base_theme + theme(legend.position = "none")

# ============================================================
# Combine with patchwork
# ============================================================

# A caption summarising replication and measurement conventions
fig_caption <- sprintf(
  paste0(
    "Each (tool, n_models, threads) cell ran in %d independent repeats; ",
    "lines and points show the median; error bars and ribbons span [min, max]. ",
    "COBRApy 2/4/8/12-thread points reflect external multiprocessing.Pool ",
    "parallelisation applied by the benchmark harness — COBRApy itself is ",
    "single-threaded per FBA call; only the 1-thread datum is its native ",
    "memory/runtime baseline (Panels C and D reflect this convention)."
  ),
  as.integer(n_repeats)
)
# Wrap the caption onto multiple lines so it never overflows the figure width
fig_caption <- paste(strwrap(fig_caption, width = 150), collapse = "\n")

# ── Panel A: fast-mic workflow schematic (horizontal, rounded boxes) ──
wf <- data.frame(
  x = 1:6,
  title = c("Inputs", "Single-species solve", "Pairwise co-culture",
            "Cross-feeding", "Classification", "Output"),
  body  = c("Genome-scale models\n(SBML)  +  medium",
            "FBA → CycleFreeFlux\n→ pFBA   (μ_alone)",
            "merge → max-min\n→ max biomass  (μ_co)",
            "donor→receiver\nexchange · GPR support",
            "benefit β, competition C\n→ 6 interaction types",
            "per-pair TSV"),
  fill = c("#E8EEF6","#DFF0E4","#DFF0E4","#FBEAD7","#F1E4F3","#ECECEC"),
  edge = c("#2166AC","#1B7837","#1B7837","#E08214","#762A83","#777777"),
  stringsAsFactors = FALSE
)
bw <- 0.90; bh <- 1.02
corners <- do.call(rbind, lapply(1:6, function(i)
  data.frame(id = i,
             x = c(i-bw/2, i+bw/2, i+bw/2, i-bw/2),
             y = c(1-bh/2, 1-bh/2, 1+bh/2, 1+bh/2),
             fill = wf$fill[i], edge = wf$edge[i], stringsAsFactors = FALSE)))
pW <- ggplot() +
  geom_segment(data = data.frame(x = 1:5),
               aes(x = x + bw/2 + 0.015, xend = x + 1 - bw/2 - 0.015, y = 1, yend = 1),
               arrow = arrow(length = unit(0.13, "cm"), type = "closed"),
               linewidth = 0.6, colour = "grey50") +
  ggforce::geom_shape(data = corners,
               aes(x = x, y = y, group = id, fill = I(fill), colour = I(edge)),
               radius = unit(1.6, "mm"), linewidth = 0.7) +
  geom_text(data = wf, aes(x = x, y = 1 + bh/2 - 0.20, label = title),
            fontface = "bold", size = 2.1, colour = "grey15") +
  geom_text(data = wf, aes(x = x, y = 1 - 0.07, label = body),
            size = 1.7, colour = "grey30", lineheight = 0.95) +
  coord_cartesian(xlim = c(0.42, 6.58), ylim = c(0.36, 1.74), expand = FALSE) +
  labs(title = "The fast-mic pipeline") +
  theme_void(base_size = 7, base_family = BMC_FONT) +
  theme(plot.title = element_text(face = "bold", size = 7.5, hjust = 0,
                                  margin = margin(b = 1)))

# Figure 1 = workflow (A) on top, then benchmark panels:
#   B = numerical agreement (accuracy, pE), C = runtime (pA),
#   D = parallel speedup (pB),              E = peak memory (pC).
# ── Panel F: fast-mic vs SMETANA — interaction-typing scalability ─────────────
# The correct interaction-tool comparator (SMETANA uses a commercial MILP solver;
# reframed has no native HiGHS). Single-thread; dashed = projection to catalogue scale.
smet_lac <- if (length(args) >= 5) args[5] else
  "results/fig1/scaling_timing_lac_L6_resistant_starch.tsv"
smet_akk <- if (length(args) >= 6) args[6] else
  "results/fig1/scaling_timing_akk_L6_resistant_starch.tsv"
SMET <- "#762A83"; N_UHGG_FULL <- 3238
humantime <- function(sec) {
  u <- c(year = 3.156e7, day = 86400, hour = 3600, min = 60)
  for (n in names(u)) if (sec >= u[[n]]) { v <- sec/u[[n]]; return(sprintf("%.0f %s%s", v, n, ifelse(v >= 2, "s", ""))) }
  sprintf("%.0f s", sec)
}
sys_of <- function(p) if (grepl("akk", tolower(basename(p)))) "Akkermansia" else "Lactobacillus"

pF <- tryCatch({
  sdat <- bind_rows(lapply(c(smet_lac, smet_akk), function(p) {
    d <- read.delim(p, stringsAsFactors = FALSE); d$system <- sys_of(p); d }))
  ssum <- sdat |> group_by(tool, system, n_pairs) |>
    summarise(wmin = min(wall_s), wmed = median(wall_s), wmax = max(wall_s), .groups = "drop") |>
    mutate(grp = paste(tool, system))
  srate <- ssum |> group_by(tool, system) |> slice_max(n_pairs, n = 1, with_ties = FALSE) |>
    ungroup() |> transmute(tool, system, rate = wmed / n_pairs, nmax = n_pairs, wnmax = wmed)
  CATALOGUE <- N_UHGG_FULL * (N_UHGG_FULL - 1) / 2
  sproj <- srate |> mutate(np = CATALOGUE, wall = rate * CATALOGUE, grp = paste(tool, system))
  sseg <- bind_rows(srate |> transmute(tool, system, np = nmax, wall = wnmax),
                    sproj |> transmute(tool, system, np, wall)) |> mutate(grp = paste(tool, system))
  refh <- data.frame(y = c(86400, 2.63e6, 3.156e7), lab = c("1 day", "1 month", "1 year"))
  tlab <- sproj |> group_by(tool) |>
    summarise(wall = max(wall), lab = sprintf("%s: %s", tool[1], humantime(max(wall))), .groups = "drop")
  ggplot(ssum, aes(n_pairs, wmed, colour = tool)) +
    geom_hline(data = refh, aes(yintercept = y), colour = GREY, linetype = "dotted",
               linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = refh, aes(x = min(ssum$n_pairs), y = y, label = lab), inherit.aes = FALSE,
              colour = GREY, size = 1.7, hjust = 0, vjust = -0.4) +
    geom_vline(xintercept = CATALOGUE, colour = GREY, linetype = "dashed", linewidth = 0.4) +
    annotate("text", x = CATALOGUE, y = min(ssum$wmed), label = sprintf("full UHGG (%s pairs)", comma(CATALOGUE)),
             colour = GREY, size = 1.6, angle = 90, vjust = -0.5, hjust = 0) +
    geom_line(data = sseg, aes(np, wall, colour = tool, linetype = system, group = grp),
              linewidth = 0.55, alpha = 0.7) +
    geom_ribbon(aes(ymin = wmin, ymax = wmax, fill = tool, group = grp), alpha = 0.13, colour = NA) +
    geom_line(aes(linetype = system, group = grp), linewidth = 0.8) +
    geom_point(size = 1) +
    geom_text(data = tlab, aes(x = CATALOGUE, y = wall, label = lab, colour = tool),
              inherit.aes = FALSE, size = 1.8, hjust = 1.03, vjust = -0.6, show.legend = FALSE) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale()),
                  expand = expansion(mult = c(0.03, 0.12))) +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "s"),
                  expand = expansion(mult = c(0.06, 0.12))) +
    annotation_logticks(sides = "bl", linewidth = 0.3) +
    scale_colour_manual(values = c("fast-mic" = BLUE, "SMETANA" = SMET), name = "Tool") +
    scale_fill_manual(values = c("fast-mic" = BLUE, "SMETANA" = SMET), name = "Tool") +
    scale_linetype_discrete(name = "System", labels = sp_md) +
    labs(x = "Number of pairs", y = "Wall time (s, log scale)",
         title = "Interaction-typing scalability",
         subtitle = "fast-mic vs SMETANA (single-thread); dashed = projection to full catalogue") +
    base_theme
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, size = 4,
                      label = "SMETANA scaling data\nnot found\n(run scaling_timing_smetana.py)") +
    theme_void(base_family = BMC_FONT) + labs(title = "Interaction-typing scalability")
})


# ── Panel G: measured peak memory for interaction typing (fast-mic vs SMETANA) ─
# Counterpart to panel E: E is the single-species benchmark (vs COBRApy), G is the
# pairwise-typing benchmark (vs SMETANA). Reads the *_mem timing runs if present.
pG <- tryCatch({
  memfiles <- c(sub("\\.tsv$", "_mem.tsv", smet_lac), sub("\\.tsv$", "_mem.tsv", smet_akk))
  memfiles <- memfiles[file.exists(memfiles)]
  stopifnot(length(memfiles) > 0)
  mdat <- bind_rows(lapply(memfiles, function(p) {
    d <- read.delim(p, stringsAsFactors = FALSE); d$system <- sys_of(p); d }))
  stopifnot("peak_rss_mb" %in% names(mdat))
  msum <- mdat |> filter(is.finite(peak_rss_mb)) |>
    group_by(tool, system, n_pairs) |>
    summarise(mb = median(peak_rss_mb), .groups = "drop")
  ratio <- msum |> tidyr::pivot_wider(names_from = tool, values_from = mb) |>
    mutate(fold = SMETANA / `fast-mic`)
  ggplot(msum, aes(factor(n_pairs), mb, fill = tool)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64,
             colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.0f", mb)), position = position_dodge(width = 0.72),
              vjust = -0.35, size = 1.55, colour = "grey15") +
    geom_text(data = ratio, aes(x = factor(n_pairs), y = SMETANA, label = sprintf("%.1f\u00d7", fold)),
              inherit.aes = FALSE, vjust = -1.8, size = 1.6, colour = GREY) +
    facet_wrap(~ system, nrow = 1, scales = "free_x", labeller = sp_labeller) +
    scale_fill_manual(values = c("fast-mic" = BLUE, "SMETANA" = SMET), name = "Tool") +
    scale_y_continuous(labels = label_number(suffix = " MB"),
                       expand = expansion(mult = c(0, 0.22))) +
    labs(x = "Number of pairs", y = "Peak resident memory",
         title = "Memory cost of interaction typing",
         subtitle = "same pairs, single process, single thread; grey = SMETANA/fast-mic ratio") +
    base_theme
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, size = 4,
                      label = "SMETANA memory data\nnot found\n(run scaling_timing_smetana.py --out-suffix _mem)") +
    theme_void(base_family = BMC_FONT) + labs(title = "Memory cost of interaction typing")
})

bench <- ((pE | pA) / (pB | pC)) + plot_layout(guides = "collect")
combined <- (pW / bench / pF / pG) +
  plot_layout(heights = c(0.75, 3, 1.5, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig <- function(path, plot, w, h, ext) {
  if (ext == "pdf") {
    ggsave(path, plot, width = w, height = h, units = "in", device = cairo_pdf)
  } else if (ext == "tiff") {
    ggsave(path, plot, width = w, height = h, units = "in", dpi = 300,
           device = "tiff", compression = "lzw")
  } else {
    ggsave(path, plot, width = w, height = h, units = "in", dpi = 300)
  }
  cat("Saved:", path, "\n")
}

for (ext in c("pdf", "png", "tiff")) {
  save_fig(paste0(out_prefix, ".", ext), combined, w = FIG_W, h = FIG_H, ext = ext)
}
