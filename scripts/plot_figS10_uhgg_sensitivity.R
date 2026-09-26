#!/usr/bin/env Rscript
# plot_figS10_uhgg_sensitivity.R  (new Supplementary Figure — answers Reviewer 1, minor: UHGG set sensitivity)
#
# The gradient pattern (inverted-U, L5 pectin peak, L6 glucose crash) is insensitive to the
# UHGG genome set: the community models already pass CheckM2 (>= 90% / <= 5%) and exclude
# Patescibacteriota (CPR); recomputing mutualism over the stricter subset that also passes these
# thresholds under the original UHGG CheckM estimates gives essentially the same curves.
#
# Input : results/figS10/uhgg_sensitivity.tsv
# Output: figS10_uhgg_sensitivity.{pdf,png,tiff}

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
INFILE <- get_arg("--in","results/figS10/uhgg_sensitivity.tsv")
OUTDIR <- get_arg("--outdir","results/figures_paper")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <- 84 / 25.4


d <- read_tsv(INFILE, show_col_types=FALSE) |>
  mutate(system=factor(system, levels=c("Akkermansia","Lactobacillus")),
         level=factor(level, levels=paste0("L",0:9))) |>
  select(system, level,
         `All 3,238 community genomes`=full_pct,
         `Also high-quality by UHGG CheckM (n = 2,903)`=strict_pct) |>
  pivot_longer(c(-system,-level), names_to="subset", values_to="mut") |>
  mutate(subset=factor(subset, levels=c("All 3,238 community genomes","Also high-quality by UHGG CheckM (n = 2,903)")))

cols <- c("All 3,238 community genomes"="#2166AC",
          "Also high-quality by UHGG CheckM (n = 2,903)"="#E6862E")

p <- ggplot(d, aes(level, mut, colour=subset, group=subset, shape=subset)) +
  geom_vline(xintercept=6.5, colour="#D6604D", linetype="dotted", linewidth=0.4) +
  geom_line(linewidth=0.5) + geom_point(size=1.3, fill="white", stroke=0.5) +
  facet_wrap(~system, scales="free_y", labeller=sp_labeller) +
  scale_colour_manual(values=cols, name=NULL) +
  scale_shape_manual(values=c(16,21), name=NULL) +
  labs(x="Prebiotic level", y="Mutualism (% of viable pairs)",
       title="The gradient pattern is insensitive to genome quality",
       subtitle="A stricter genome-quality filter preserves the L5 peak and the L6 glucose crash (dotted red = L5→L6)") +
  theme_classic(base_size=8, base_family = BMC_FONT) +
  theme(legend.position="top", legend.text=element_text(size=6.5),
        legend.key.size=unit(0.32,"cm"),
        axis.title=element_text(size=7), axis.text=element_text(size=6),
        plot.title=element_text(face="bold", size=8.5),
        plot.subtitle=element_text(size=5.6, colour="grey30"),
        strip.text=md(face="bold", size=7))

for(e in c("pdf","png","tiff")){
  fp<-file.path(OUTDIR,paste0("figS10_uhgg_sensitivity.",e))
  if(e=="tiff") ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300,device="tiff",compression="lzw")
  else if(e=="pdf") ggsave(fp,p,width=FIG_W,height=FIG_H, device=cairo_pdf) else ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300)
}
cat("Saved figS10_uhgg_sensitivity to",OUTDIR,"\n")
