## figure_theme.R — shared figure conventions for the fast-mic manuscript
##
## BMC / Microbiome figure requirements (submission guidelines):
##   * image resolution ~300 dpi at final size
##   * full-page width 170 mm, half-page width 85 mm, maximum height 225 mm
##     (figure plus legend)
##   * every line wider than 0.25 pt
##   * all text legible at the final size
## BMC does not prescribe a font family; we use Arial throughout for consistency
## and fall back to the generic sans family when Arial is unavailable.
##
## Species names must be italic (ASM/ICNP convention, and BMC copy-editing
## applies it in the text). ggplot draws labels as plain strings, so the helpers
## below convert species names to markdown and the theme renders that markdown
## with ggtext. Data values are never modified — only how they are drawn.
##
## Usage in a plotting script:
##   source("scripts/figure_theme.R")
##   ... + labs(title = sp_md("Akkermansia cooperates ...")) +
##       scale_x_discrete(labels = sp_md) +
##       facet_wrap(~ system, labeller = sp_labeller) +
##       theme_bmc()            # or: your_theme + theme(<element> = element_markdown())

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
})

## ── page geometry (mm -> inches), for ggsave calls ───────────────────────────
BMC_FULL_WIDTH_IN <- 170 / 25.4   # 6.69 in
BMC_HALF_WIDTH_IN <- 85 / 25.4    # 3.35 in
BMC_MAX_HEIGHT_IN <- 225 / 25.4   # 8.86 in
BMC_DPI <- 300
BMC_MIN_LINEWIDTH <- 0.25 / .pt   # 0.25 pt expressed in ggplot linewidth units

## ── font family, with a graceful fallback ────────────────────────────────────
.bmc_font <- local({
  fam <- "Arial"
  ok <- tryCatch({
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      any(systemfonts::system_fonts()$family == fam)
    } else TRUE
  }, error = function(e) FALSE)
  if (ok) fam else ""      # "" = device default sans
})
BMC_FONT <- .bmc_font

## Text drawn by geoms (geom_text, geom_label, ggrepel) does not inherit the
## theme font, so set it once here for every script that sources this file.
local({
  for (g in c("text", "label", "text_repel", "label_repel")) {
    try(ggplot2::update_geom_defaults(g, list(family = BMC_FONT)), silent = TRUE)
  }
})

## ── species names ────────────────────────────────────────────────────────────
## Genus names that occur in this study, plus the genera named in figure
## annotations. Epithets are italicised together with the genus.
.SPECIES_GENERA <- c(
  "Akkermansia", "Lactobacillus", "Lacticaseibacillus", "Limosilactobacillus",
  "Lactiplantibacillus", "Ligilactobacillus", "Leuconostoc", "Lactococcus",
  "Bifidobacterium", "Eubacterium", "Anaerobutyricum", "Anaerostipes",
  "Faecalibacterium", "Bacteroides", "Ruminococcus", "Blautia", "Lachnospira",
  "Phascolarctobacterium", "Salmonella", "Escherichia", "Clostridium",
  "Roseburia", "Subdoligranulum", "Latilactobacillus"
)

## Species epithets that occur in this study (plus the gut taxa named in figure
## annotations). A genus is only joined to the following word when that word is
## on this list — otherwise "Lactobacillus grows" and "Akkermansia and" would
## italicise an ordinary English word along with the genus.
.SPECIES_EPITHETS <- c(
  "muciniphila", "massiliensis", "biwaensis", "glycaniphila",
  "acidophilus", "gasseri", "crispatus", "delbrueckii", "casei", "paracasei",
  "rhamnosus", "plantarum", "reuteri", "johnsonii", "helveticus", "salivarius",
  "fermentum", "brevis", "longum", "breve", "adolescentis", "bifidum",
  "thetaiotaomicron", "fragilis", "ovatus", "vulgatus", "uniformis",
  "xylanisolvens", "prausnitzii", "rectale", "hallii", "gnavus", "bromii",
  "coli", "difficile", "enterica", "intestinalis"
  ## "sp." / "spp." are abbreviations, not epithets: they stay roman
)

#' Italicise species names inside a label, returning markdown.
#'
#' Handles three forms:
#'   "Akkermansia"              -> "*Akkermansia*"
#'   "Akkermansia muciniphila"  -> "*Akkermansia muciniphila*"
#'   "L. gasseri ATCC33323"     -> "*L. gasseri* ATCC33323"   (strain code stays roman)
#' A genus followed by an ordinary word keeps that word roman:
#'   "Lactobacillus grows"      -> "*Lactobacillus* grows"
#' Anything already wrapped in asterisks is left alone, and non-character input
#' is returned unchanged so the helpers are safe to pass to scale labels.
sp_md <- function(x) {
  if (!is.character(x)) return(x)
  out <- x
  genus   <- paste0("\\b(", paste(.SPECIES_GENERA, collapse = "|"), ")\\b")
  epithet <- paste0("(\\s+(?:", paste(.SPECIES_EPITHETS, collapse = "|"), ")\\b\\.?)?")
  ## genus (+ following epithet, only when it really is one)
  out <- gsub(paste0(genus, epithet), "*\\1\\2*", out, perl = TRUE)
  ## abbreviated genus + epithet, e.g. "A. muciniphila", "L. gasseri"
  out <- gsub("(?<![*\\w])([A-Z]\\.)\\s([a-z]{3,})\\b", "*\\1 \\2*", out, perl = TRUE)
  ## tidy up doubled markers produced by adjacent matches
  out <- gsub("\\*{2,}", "*", out)
  out
}

#' Italicise the species part of a label as a plotmath expression string.
#'
#' For geoms that cannot render markdown (ggtree::geom_tiplab and anything else
#' drawn with parse = TRUE). Only the leading species name is italicised; a
#' trailing strain designation or accession stays roman, as the nomenclature
#' rules require:
#'   "A. muciniphila GCF_000020225.1" -> italic("A. muciniphila")*" GCF_000020225.1"
sp_plotmath <- function(x) {
  gen <- paste0("(", paste(.SPECIES_GENERA, collapse = "|"), "|[A-Z]\\.)")
  pat <- paste0("^", gen, "(\\s+[a-z]{3,})?")
  vapply(as.character(x), function(s) {
    m <- regexpr(pat, s, perl = TRUE)
    if (m == -1) return(paste0('"', s, '"'))
    len <- attr(m, "match.length")
    sp   <- substr(s, 1, len)
    rest <- substring(s, len + 1)
    if (nzchar(trimws(rest))) sprintf('italic("%s")*"%s"', sp, rest)
    else                      sprintf('italic("%s")', sp)
  }, character(1), USE.NAMES = FALSE)
}

#' facet labeller that italicises species names
sp_labeller <- ggplot2::as_labeller(sp_md)

#' Convenience: element_markdown() carrying the BMC font
md <- function(...) ggtext::element_markdown(family = BMC_FONT, ...)

## ── theme ────────────────────────────────────────────────────────────────────
#' BMC-compliant base theme. Text elements that commonly carry species names
#' (title, subtitle, axis text, strip text, legend text) are rendered as
#' markdown so that sp_md() output appears in italics.
theme_bmc <- function(base_size = 8, base_family = BMC_FONT, grid = FALSE) {
  base <- if (grid) ggplot2::theme_bw(base_size = base_size, base_family = base_family)
          else       ggplot2::theme_classic(base_size = base_size, base_family = base_family)
  base + ggplot2::theme(
    plot.title      = md(size = base_size + 1.5, face = "bold"),
    plot.subtitle   = md(size = base_size - 0.5, colour = "grey30"),
    axis.text       = md(size = base_size - 0.5),
    axis.title      = ggplot2::element_text(size = base_size, family = base_family),
    strip.text      = md(size = base_size, face = "bold"),
    legend.text     = md(size = base_size - 0.5),
    legend.title    = ggplot2::element_text(size = base_size - 0.5, family = base_family),
    plot.tag        = ggplot2::element_text(size = base_size + 2, face = "bold", family = base_family),
    line            = ggplot2::element_line(linewidth = max(BMC_MIN_LINEWIDTH, 0.3))
  )
}

#' Add markdown rendering (and the BMC font) to an existing plot without
#' replacing its theme — used where a script already has a tuned theme.
add_md <- function(p, title = TRUE, subtitle = TRUE, axis = TRUE,
                   strip = TRUE, legend = TRUE) {
  th <- list()
  if (title)    th$plot.title    <- md(face = "bold")
  if (subtitle) th$plot.subtitle <- md(colour = "grey30")
  if (axis)     th$axis.text     <- md()
  if (strip)    th$strip.text    <- md()
  if (legend)   th$legend.text   <- md()
  p + do.call(ggplot2::theme, th)
}

#' ggsave wrapper that enforces the BMC geometry and resolution.
ggsave_bmc <- function(filename, plot, width_mm = 170, height_mm = 225, ...) {
  stopifnot(width_mm <= 170, height_mm <= 225)
  ext <- tolower(tools::file_ext(filename))
  args <- list(filename = filename, plot = plot,
               width = width_mm / 25.4, height = height_mm / 25.4,
               units = "in", dpi = BMC_DPI, ...)
  if (ext == "pdf") args$device <- grDevices::cairo_pdf
  if (ext %in% c("tif", "tiff")) { args$device <- "tiff"; args$compression <- "lzw" }
  do.call(ggplot2::ggsave, args)
}
