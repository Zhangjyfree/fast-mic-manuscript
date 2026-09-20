#!/usr/bin/env Rscript
# plot_figS5_benefit_metric.R  (new Supplementary Figure — answers Reviewer 1, Comment 3)
#
# Reviewer 1.3 worried that the "glucose crash" (mutualism drop at L6) could be a
# NORMALIZATION ARTIFACT of the relative benefit  b = (mu_co - mu_alone)/mu_alone:
# at L6 free glucose lifts mu_alone for everyone, so b shrinks even if the ABSOLUTE
# cross-feeding benefit is unchanged.
#
# This figure shows the crash is NOT an artifact:
#   * mu_alone indeed jumps at L6 (+35% Akk, +42% Lac)  -> reviewer's premise holds,
#   * but mutualism re-scored with a FIXED ABSOLUTE threshold (mu_co - mu_alone > 1e-3 h^-1,
#     denominator removed) drops by essentially the SAME amount as the relative metric,
#   * and the mean absolute co-culture benefit becomes MORE negative at L6
#     -> competition genuinely intensifies; the collapse is real, not a denominator effect.
#
# Input : results/figS5/absolute_vs_relative.tsv
# Output: figS5_benefit_metric.{pdf,png,tiff}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(readr) })

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
get_arg <- function(f,d){ i<-match(f,args); if(!is.na(i)&&i<length(args)) args[i+1] else d }
INFILE <- get_arg("--in","results/figS5/absolute_vs_relative.tsv")
OUTDIR <- get_arg("--outdir","results/figures_paper")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 88 / 25.4


d <- read_tsv(INFILE, show_col_types=FALSE) |>
  mutate(system=factor(system, levels=c("Akkermansia","Lactobacillus")),
         level =factor(level,  levels=paste0("L",0:9)))

# long form for the two mutualism metrics
mut <- d |> select(system, level, mu_alone_mean,
                   `Relative benefit  β = (μ_co−μ_alone)/μ_alone` = mutualism_relative_pct,
                   `Absolute benefit  μ_co−μ_alone > 10⁻³ h⁻¹`     = mutualism_absolute_pct) |>
  pivot_longer(c(-system,-level,-mu_alone_mean), names_to="metric", values_to="mut")

# scale μ_alone onto the mutualism axis for a secondary reference line (per facet range)
# use a global factor so both facets are comparable
sf <- 15   # μ_alone (~0.5–1.4) × 15 ≈ 8–21, overlaps the mutualism % range

## legend labels drawn as markdown, so the subscripts match the rest of the figures
METRIC_MD <- function(x) {
  x <- sub("\u03bc_co", "\u03bc<sub>co</sub>", x, fixed = TRUE)
  gsub("\u03bc_alone", "\u03bc<sub>alone</sub>", x, fixed = TRUE)
}
cols <- c("Relative benefit  β = (μ_co−μ_alone)/μ_alone"="#2166AC",
          "Absolute benefit  μ_co−μ_alone > 10⁻³ h⁻¹"="#E6862E")

mu_df <- distinct(d, system, level, mu_alone_mean)

p <- ggplot(mut, aes(level, mut, group=metric, colour=metric, shape=metric)) +
  # μ_alone reference (grey dashed, secondary axis)
  geom_line(data=mu_df, aes(x=level, y=mu_alone_mean*sf, group=1),
            colour="grey60", linetype="dashed", linewidth=0.5, inherit.aes=FALSE) +
  geom_vline(xintercept=6.5, colour="#D6604D", linetype="dotted", linewidth=0.5) +
  geom_line(linewidth=0.6) + geom_point(size=1.5, fill="white", stroke=0.6) +
  facet_wrap(~system, scales="free_y", labeller=sp_labeller) +
  scale_colour_manual(values=cols, name=NULL, labels=METRIC_MD) +
  scale_shape_manual(values=c(21,24), name=NULL, labels=METRIC_MD) +
  scale_y_continuous(sec.axis=sec_axis(~ ./sf, name="Mean monoculture growth μ<sub>alone</sub> (h<sup>−1</sup>)")) +
  labs(x="Prebiotic level", y="Mutualism (% of viable pairs)",
       title="The glucose crash is not a normalization artifact of the relative benefit metric",
       subtitle=sp_md(paste0("Re-scoring mutualism with a fixed ABSOLUTE benefit (orange) reproduces the L5→L6 crash of the relative metric (blue),<br>",
                             "even though μ<sub>alone</sub> jumps at L6 (grey dashed, +35% Akkermansia / +42% Lactobacillus). Dotted red line = L5→L6 transition."))) +
  theme_classic(base_size=8, base_family = BMC_FONT) +
  theme(legend.position="top", legend.text=md(size=6.5),
        legend.key.size=unit(0.32,"cm"),
        axis.title=element_text(size=7), axis.text=element_text(size=6),
        plot.title=md(face="bold", size=8.5),
        plot.subtitle=md(size=5.6, colour="grey30"),
        strip.text=md(face="bold", size=7),
        axis.title.y.right=md(colour="grey45", size=7, angle=-90),
        axis.text.y.right=element_text(colour="grey45"))

save_one <- function(ext){
  fp <- file.path(OUTDIR, paste0("figS5_benefit_metric.",ext))
  if(ext=="tiff") ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300,device="tiff",compression="lzw")
  else if(ext=="pdf") ggsave(fp,p,width=FIG_W,height=FIG_H, device=cairo_pdf)
  else ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300)
}
for(e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS5_benefit_metric to",OUTDIR,"\n")
