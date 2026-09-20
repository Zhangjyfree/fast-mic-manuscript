#!/usr/bin/env Rscript
# plot_figS4_strain_bootstrap.R  (new Supplementary Figure — answers Reviewer 1, Comment 6)
#
# Strain-level statistics for the two headline contrasts, with the probiotic STRAIN
# as the unit of replication (bootstrap across strains; Akkermansia n=6, Lactobacillus n=10):
#   (A) Generalist vs specialist at L5 — per-strain mutualism %, genus mean ± 95% bootstrap CI.
#       Lactobacillus/Akkermansia ratio = 2.01× (95% CI 1.34–3.12); difference 10.9 pts
#       all statistics (ratio, difference, CIs, p-values) come from strain_level_tests.tsv.
#   (B) L5→L6 collapse — paired per-strain change; every strain declines.
#       Akkermansia −3.7 pts (95% CI −5.2 to −2.2), Wilcoxon p = 0.031;
#       Lactobacillus −5.1 pts (95% CI −7.5 to −2.6), Wilcoxon p = 0.0039.
#
# Input : results/figS4/per_strain_mutualism_L5_L6.tsv
#         results/figS4/genus_ci_summary.tsv
#         results/figS4/strain_level_tests.tsv      <- all statistics, from stats_strain_level.py
#         (run `python3 scripts/stats_strain_level.py` to regenerate all three)
# Output: figS4_strain_bootstrap.{pdf,png,tiff}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork) })

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
DIR    <- get_arg("--dir","results/figS4")
OUTDIR <- get_arg("--outdir","results/figures_paper")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 82 / 25.4


ps  <- read_tsv(file.path(DIR,"per_strain_mutualism_L5_L6.tsv"), show_col_types=FALSE) |>
  mutate(system=recode(system, akk="Akkermansia", lac="Lactobacillus"))
st  <- read_tsv(file.path(DIR,"strain_level_tests.tsv"), col_types=cols(.default="c"))
getst <- function(pat, col) { v <- st[[col]][grepl(pat, st$quantity, fixed=TRUE)][1]; v }
lab_A <- sprintf("Lac/Akk = %s\u00d7 (95%% CI %s\u2013%s)\nMann\u2013Whitney p = %s",
                 getst("mutualism ratio at L5","estimate"),
                 getst("mutualism ratio at L5","ci_lo"), getst("mutualism ratio at L5","ci_hi"),
                 getst("rank test across strains","p_value"))
ci  <- read_tsv(file.path(DIR,"genus_ci_summary.tsv"), show_col_types=FALSE) |>
  mutate(system=factor(system, levels=c("Akkermansia","Lactobacillus")))

cols <- c("Akkermansia"="#762A83","Lactobacillus"="#1B7837")

# ---- Panel A: L5 genus comparison ----
ciA <- ci |> filter(level=="L5") |> mutate(system=factor(system,levels=c("Akkermansia","Lactobacillus")))
psA <- ps |> mutate(system=factor(system,levels=c("Akkermansia","Lactobacillus")))
pA <- ggplot(ciA, aes(system, mean, colour=system)) +
  geom_jitter(data=psA, aes(system, L5_mut_pct, colour=system), width=0.12, size=1.4, alpha=0.55) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.16, linewidth=0.8) +
  geom_point(size=2.2) +
  annotate("text", x=1.5, y=36, label=lab_A, size=2, colour="grey25", lineheight=0.95) +
  scale_colour_manual(values=cols, guide="none") +
  scale_x_discrete(labels=sp_md) +
  scale_y_continuous(limits=c(0,38), expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL, y="Mutualism at L5 (% of viable pairs)",
       title="A. Generalist vs specialist (strain-level, L5)",
       subtitle="Each point = one probiotic strain; bar = genus mean ± 95% bootstrap CI") +
  theme_classic(base_size=8, base_family = BMC_FONT) +
  theme(plot.title=element_text(face="bold",size=8.5),
        axis.title=element_text(size=7), axis.text=element_text(size=6),
        plot.subtitle=element_text(size=5.6,colour="grey30"),
        axis.text.x=md(face="bold",size=6.5))

# ---- Panel B: L5 -> L6 paired change ----
psL <- ps |> pivot_longer(c(L5_mut_pct,L6_mut_pct), names_to="level", values_to="mut") |>
  mutate(level=recode(level, L5_mut_pct="L5", L6_mut_pct="L6"),
         system=factor(system,levels=c("Akkermansia","Lactobacillus")))
ciB <- ci |> mutate(level=factor(level,levels=c("L5","L6")))
lab <- tibble(system=factor(c("Akkermansia","Lactobacillus"),levels=c("Akkermansia","Lactobacillus")),
              txt=c("−3.7 pts\nWilcoxon p = 0.031","−5.1 pts\nWilcoxon p = 0.0039"), x=1.5, y=c(2,4))
pB <- ggplot(psL, aes(level, mut, group=strain, colour=system)) +
  geom_line(alpha=0.35, linewidth=0.5) + geom_point(alpha=0.5, size=1.1) +
  geom_line(data=ciB, aes(level, mean, group=1), colour="grey20", linewidth=1.1, inherit.aes=FALSE) +
  geom_errorbar(data=ciB, aes(level, ymin=lo, ymax=hi), width=0.1, colour="grey20", linewidth=0.7, inherit.aes=FALSE) +
  geom_point(data=ciB, aes(level, mean), colour="grey20", size=1.8, inherit.aes=FALSE) +
  geom_text(data=lab, aes(x=x, y=y, label=txt), inherit.aes=FALSE, size=2, colour="grey25", lineheight=0.95) +
  facet_wrap(~system, scales="free_y", labeller=sp_labeller) +
  scale_colour_manual(values=cols, guide="none") +
  labs(x=NULL, y="Mutualism (% of viable pairs)",
       title="B. L5→L6 collapse is paired and per-strain",
       subtitle=sp_md("Thin lines = individual strains (all 6 Akkermansia and 9/10 Lactobacillus decline);<br>thick grey = genus mean ± 95% bootstrap CI")) +
  theme_classic(base_size=8, base_family = BMC_FONT) +
  theme(plot.title=element_text(face="bold",size=8.5),
        axis.title=element_text(size=7), axis.text=element_text(size=6),
        plot.subtitle=md(size=5.6,colour="grey30"),
        strip.text=md(face="bold",size=7))

p <- pA | pB
save_one <- function(ext){
  fp <- file.path(OUTDIR, paste0("figS4_strain_bootstrap.",ext))
  if(ext=="tiff") ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300,device="tiff",compression="lzw")
  else if(ext=="pdf") ggsave(fp,p,width=FIG_W,height=FIG_H, device=cairo_pdf)
  else ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300)
}
for(e in c("pdf","png","tiff")) save_one(e)
cat("Saved figS4_strain_bootstrap to",OUTDIR,"\n")
