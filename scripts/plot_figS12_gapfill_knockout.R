#!/usr/bin/env Rscript
# plot_figS12_gapfill_knockout.R  (Supplementary Figure S12 — answers Reviewer 1.5 / 2.5)
# Predicted cooperation is robust to removal of gap-fill-candidate reactions
# (direct classification-level test).
#
# Input : results/figS12/R15_summary.tsv
# Output: results/figures_paper/figS12_gapfill_knockout.{pdf,png,tiff}
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
INFILE <- get_arg("--in","results/figS12/R15_summary.tsv")
OUTDIR <- get_arg("--outdir","results/figures_paper")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

FIG_W <- 170 / 25.4      # BMC full-page width
FIG_H <-  85 / 25.4

d <- read_tsv(INFILE, show_col_types=FALSE) |>
  mutate(system=factor(system, levels=c("Akkermansia","Lactobacillus")))
cols <- c("Akkermansia"="#762A83","Lactobacillus"="#1B7837")

# Panel A: overall mutualism % before vs after
pa_df <- d |> select(system, Control=ctrl_mut_pct, `Gap-fill KO`=aprime_mut_pct) |>
  pivot_longer(-system, names_to="condition", values_to="mut") |>
  mutate(condition=factor(condition, levels=c("Control","Gap-fill KO")))
pA <- ggplot(pa_df, aes(condition, mut, fill=system)) +
  geom_col(width=0.6, colour="grey25", linewidth=0.25) +
  geom_text(aes(label=sprintf("%.1f%%",mut)), vjust=-0.4, size=2.2) +
  facet_wrap(~system, scales="free_y", labeller=sp_labeller) +
  scale_fill_manual(values=cols, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL, y="Mutualism (% of viable pairs)",
       title="A. Overall mutualism is nearly unchanged",
       subtitle="After blocking dormant no-GPR internal reactions\n(transporters kept; monoculture growth unchanged)") +
  theme_classic(base_size=8, base_family = BMC_FONT) +
  theme(plot.title=element_text(face="bold",size=8.5), plot.subtitle=element_text(size=6,colour="grey30"),
        axis.title=element_text(size=7), axis.text=element_text(size=6),
        strip.text=md(face="bold",size=7))

# Panel B: fate of control-mutualistic pairs (stacked)
pb_df <- d |> transmute(system,
                        `Retained mutualism`=retained_mut,
                        `→ competition`=to_competition,
                        `→ neutral / commensal`=to_neutral_commensal) |>
  pivot_longer(-system, names_to="fate", values_to="n") |>
  group_by(system) |> mutate(pct=100*n/sum(n)) |> ungroup() |>
  mutate(fate=factor(fate, levels=c("→ neutral / commensal","→ competition","Retained mutualism")))
fcols <- c("Retained mutualism"="#4DAC26","→ competition"="#D6604D","→ neutral / commensal"="#BBBBBB")
pB <- ggplot(pb_df, aes(system, pct, fill=fate)) +
  geom_col(width=0.6, colour="grey25", linewidth=0.25) +
  geom_text(data=pb_df |> filter(fate=="Retained mutualism"),
            aes(label=sprintf("%.1f%%",pct)), position=position_stack(vjust=0.5), size=2.4, fontface="bold", colour="white") +
  scale_fill_manual(values=fcols, name=NULL) +
  scale_x_discrete(labels=sp_md) +
  scale_y_continuous(expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL, y="Fate of control-mutualistic pairs (%)",
       title="B. ≥90% of mutualistic pairs keep their classification",
       subtitle=sp_md("Akkermansia n=1,655; Lactobacillus n=5,753<br>control-mutualistic viable pairs")) +
  theme_classic(base_size=8, base_family = BMC_FONT) +
  theme(legend.position="top", legend.text=element_text(size=6.5),
        legend.key.size=unit(0.32,"cm"),
        plot.title=element_text(face="bold",size=8.5), plot.subtitle=md(size=6,colour="grey30"),
        axis.title=element_text(size=7), axis.text=element_text(size=6),
        axis.text.x=md(face="bold",size=6.5))

p <- pA | pB
for(e in c("pdf","png","tiff")){
  fp<-file.path(OUTDIR,paste0("figS12_gapfill_knockout.",e))
  if(e=="tiff") ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300,device="tiff",compression="lzw")
  else if(e=="pdf") ggsave(fp,p,width=FIG_W,height=FIG_H, device=cairo_pdf) else ggsave(fp,p,width=FIG_W,height=FIG_H,dpi=300)
}
cat("Saved figS12_gapfill_knockout (PDF/PNG/TIFF) to", OUTDIR, "\n")
