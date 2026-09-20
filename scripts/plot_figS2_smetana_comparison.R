#!/usr/bin/env Rscript
# Fig S13: fast-mic vs SMETANA — projected cost + concordance (defensive supplement).
#   A  Projected wall time at full-catalogue scale (bars; complements main Fig 1F)
#   B  Do SMETANA's set-based POTENTIAL scores (MIP/MRO) rank fast-mic's REALIZED
#      calls? MIP / MRO by fast-mic interaction type, with rank AUROC.
#   C  Rank scatters: fast-mic cross-feeding # vs MIP, competition C vs MRO (Spearman).
# Honest framing: MIP/MRO are set-based potential, fast-mic C/cross-feeding are flux-based
# realized metrics, so tight agreement is not expected; AUROC~0.5 / rho~0 means the two
# tools capture complementary aspects, not that either is wrong.
#
# Inputs: results/fig1/scaling_timing_{lac,akk}_L6_resistant_starch.tsv (A)
#         results/figS2/smetana_{lac,akk}_{L5_pectin,L6_resistant_starch}.tsv (B,C)
# Usage: Rscript scripts/plot_figS2_smetana_comparison.R [output_prefix]

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(scales)
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

FM <- "."
args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[1] else "results/figures_paper/figS2_smetana_comparison"
N_UHGG_FULL <- 3238; N_MEDIA <- 10
BLUE <- "#2166AC"; SMET <- "#762A83"; GREY <- "#999999"; MUT <- "#1B7837"; COMP <- "#B2182B"
tool_col <- c("fast-mic" = BLUE, "SMETANA" = SMET); type_col <- c("mutualism" = MUT, "competition" = COMP)

base_theme <- theme_bw(base_size = 7.5, base_family = BMC_FONT) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5), axis.text = element_text(size = 6),
        axis.title = element_text(size = 7), plot.title = element_text(size = 8, face = "bold"),
        strip.text = md(size = 7, face = "bold"),
        plot.subtitle = element_text(size = 5.8, colour = "grey40"))

humantime <- function(sec) { u <- c(year=3.156e7, day=86400, hour=3600, min=60)
  for (n in names(u)) if (sec >= u[[n]]) { v <- sec/u[[n]]; return(sprintf("%.0f %s%s", v, n, ifelse(v>=2,"s",""))) }
  sprintf("%.0f s", sec) }
auroc <- function(score, pos) { ok <- is.finite(score) & !is.na(pos); score<-score[ok]; pos<-pos[ok]
  n1<-sum(pos); n2<-sum(!pos); if(n1==0||n2==0) return(NA_real_)
  r<-rank(score); (sum(r[pos]) - n1*(n1+1)/2)/(n1*n2) }
sp_rho <- function(x,y){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<3) return(NA_real_)
  suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
sys_of <- function(p) if (grepl("akk", tolower(basename(p)))) "Akkermansia" else "Lactobacillus"

# ── Panel A: projected-cost bars (from scaling timing) ───────────────────────
scal <- bind_rows(lapply(
  c(sprintf("%s/results/fig1/scaling_timing_lac_L6_resistant_starch.tsv", FM),
    sprintf("%s/results/fig1/scaling_timing_akk_L6_resistant_starch.tsv", FM)),
  function(p) { d <- read.delim(p, stringsAsFactors=FALSE); d$system <- sys_of(p); d }))
rate <- scal |> group_by(tool, system, n_pairs) |> summarise(wmed=median(wall_s), .groups="drop") |>
  group_by(tool, system) |> slice_max(n_pairs, n=1, with_ties=FALSE) |> ungroup() |>
  mutate(rate = wmed/n_pairs, n_prob = ifelse(system=="Akkermansia", 6, 10))
CATALOGUE <- N_UHGG_FULL*(N_UHGG_FULL-1)/2
barp <- rate |> tidyr::crossing(scale = c("This study","Full UHGG\nall-vs-all")) |>
  mutate(npair = ifelse(scale=="This study", n_prob*N_UHGG_FULL*N_MEDIA, CATALOGUE),
         wall = rate*npair, lab = vapply(wall, humantime, character(1)))
pA <- ggplot(barp, aes(scale, wall, fill=tool)) +
  geom_col(position=position_dodge(width=0.75), width=0.68) +
  geom_text(aes(label=lab), position=position_dodge(width=0.75), vjust=-0.4, size=1.8, colour="grey15") +
  facet_wrap(~system, nrow=1, labeller=sp_labeller) +
  scale_y_log10(labels=label_number(scale_cut=cut_short_scale(), suffix="s"), expand=expansion(mult=c(0,0.22))) +
  scale_fill_manual(values=tool_col, name="Tool") +
  labs(x=NULL, y="Projected wall time (s, log)", title="Projected cost at catalogue scale",
       subtitle="median s/pair × pair count (single-thread)") + base_theme


# ── B/C: concordance (read whatever figS2 batch tsvs exist) ──────────────────
cfiles <- list.files(sprintf("%s/results/figS2", FM), pattern="^smetana_.*\\.tsv$", full.names=TRUE)
num <- function(v) suppressWarnings(as.numeric(ifelse(v %in% c("None","NA",""), NA, v)))
build_concord <- function() {
  d <- bind_rows(lapply(cfiles, function(p){
    x <- read_tsv(p, show_col_types=FALSE)
    x$system <- sys_of(p); x$level <- if(grepl("L6",p)) "L6" else "L5"; x }))
  d <- d |> mutate(smetana_MIP=num(smetana_MIP), smetana_MRO=num(smetana_MRO),
                   fm_C=num(fm_C), fm_xfeed=num(fm_xfeed)) |>
    filter(fm_type %in% c("mutualism","competition"))
  is_mut <- d$fm_type=="mutualism"
  a_mip <- auroc(d$smetana_MIP, is_mut); a_mro <- auroc(d$smetana_MRO, !is_mut)
  r_mip <- sp_rho(d$fm_xfeed, d$smetana_MIP); r_mro <- sp_rho(d$fm_C, d$smetana_MRO)
  labB <- data.frame(var=c("SMETANA MIP\n(cooperation potential)","SMETANA MRO\n(competition potential)"),
                     lab=c(sprintf("AUROC=%.2f", a_mip), sprintf("AUROC=%.2f", a_mro)))
  dB <- bind_rows(
    d |> transmute(fm_type, value=smetana_MIP, var="SMETANA MIP\n(cooperation potential)"),
    d |> transmute(fm_type, value=smetana_MRO, var="SMETANA MRO\n(competition potential)"))
  pB <- ggplot(dB, aes(fm_type, value, fill=fm_type)) +
    geom_violin(alpha=0.35, colour=NA, scale="width") + geom_boxplot(width=0.22, outlier.size=0.3, alpha=0.9) +
    facet_wrap(~var, scales="free_y") +
    geom_text(data=labB, aes(x=1.5, y=Inf, label=lab), inherit.aes=FALSE, vjust=1.3, size=2, colour="grey20") +
    scale_fill_manual(values=type_col, guide="none") +
    labs(x=NULL, y="SMETANA potential score", title="Potential vs realized: does SMETANA rank fast-mic's calls?",
         subtitle=sprintf("n=%d pairs (lac+akk, L5+L6); mutualism should sit higher on MIP, competition on MRO", nrow(d))) +
    base_theme
  pC1 <- ggplot(d, aes(fm_xfeed, smetana_MIP, colour=fm_type)) + geom_point(alpha=0.45, size=0.7) +
    annotate("text", x=-Inf, y=Inf, label=sprintf("rho=%.2f", r_mip), hjust=-0.15, vjust=1.5, size=2, colour="grey20") +
    scale_colour_manual(values=type_col, name=NULL) +
    labs(x="fast-mic cross-fed metabolites (#)", y="SMETANA MIP", title="Realized cross-feeding vs MIP") + base_theme
  pC2 <- ggplot(d, aes(fm_C, smetana_MRO, colour=fm_type)) + geom_point(alpha=0.45, size=0.7) +
    annotate("text", x=-Inf, y=Inf, label=sprintf("rho=%.2f", r_mro), hjust=-0.15, vjust=1.5, size=2, colour="grey20") +
    scale_colour_manual(values=type_col, name=NULL) +
    labs(x="fast-mic competition intensity C", y="SMETANA MRO", title="Realized competition vs MRO") + base_theme
  list(pB=pB, pC=(pC1 | pC2))
}

if (length(cfiles) > 0) {
  bc <- build_concord()
  combined <- (pA / bc$pB / bc$pC) + plot_layout(heights=c(1, 1.1, 1)) +
    plot_annotation(tag_levels=list(c("A","B","C","")))   # C tags the first plot of the row
} else {
  msg <- ggplot() + annotate("text", x=.5, y=.5, size=4,
    label="Concordance data not found\n(run run_smetana_batch.py → results/figS2/)") + theme_void(base_family = BMC_FONT)
  combined <- (pA / msg) + plot_layout(heights=c(1,1.2)) + plot_annotation(tag_levels="A")
}

for (ext in c("pdf","png","tiff")) {
  path <- paste0(out_prefix, ".", ext)
  FIG_W <- 170 / 25.4      # BMC full-page width
  FIG_H <- 195 / 25.4      # within the 225 mm limit
  if (ext=="tiff") ggsave(path, combined, width=FIG_W, height=FIG_H, units="in", dpi=300, device="tiff", compression="lzw")
  else if (ext=="png") ggsave(path, combined, width=FIG_W, height=FIG_H, units="in", dpi=300)
  else ggsave(path, combined, width=FIG_W, height=FIG_H, units="in", device=cairo_pdf)
  cat("Saved:", path, "\n")
}
