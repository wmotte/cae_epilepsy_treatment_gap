# ==============================================================================
# plot_theme.R
# Shared figure house style for the Lancet Neurology submission.
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Defines exactly:
#   bc_palette    : named Okabe-Ito colourblind-safe palette
#   state_colours : stable health-state -> colour map (reused by every state fig)
#   country_colours : stable country -> colour map (Nigeria/Ecuador/UK)
#   theme_cae()   : theme_classic base + top legend + bold facet strips
#   number_ticks(): pretty()-based axis-tick helper
#   save_fig()    : writes BOTH a 600-dpi PNG and a vector PDF (cairo_pdf)
#
# Follows the repository-level FIGURE_STYLE.md.
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

# ---- Okabe-Ito palette (colourblind-safe). Preference order for series:
#   orange (1st) -> skyblue (2nd) -> green (3rd) -> blue (4th) -> vermillion ...
bc_palette <- c(grey = "#999999", orange = "#E69F00", skyblue = "#56B4E9",
                green = "#009E73", blue = "#0072B2", vermillion = "#D55E00",
                purple = "#CC79A7", yellow = "#F0E442")

# ---- Stable health-state colour mapping (used by fig1, fig8, figS1, figS3,
# figS4 - every figure that encodes states). Ordered sickest -> healthiest so
# the palette reads as a clinical gradient, with grey reserved for Death.
state_colours <- c(
  Untreated  = unname(bc_palette["vermillion"]),  # worst: symptomatic, no ASM
  Treated    = unname(bc_palette["orange"]),       # on ASM, still seizing
  Controlled = unname(bc_palette["skyblue"]),      # on ASM, seizure-free
  Remission  = unname(bc_palette["green"]),         # off ASM, seizure-free
  Death      = unname(bc_palette["grey"])           # absorbing
)

# ---- Stable country colour mapping (three series -> 1st/2nd/3rd preference).
country_colours <- c(
  Nigeria = unname(bc_palette["orange"]),
  Ecuador = unname(bc_palette["skyblue"]),
  UK      = unname(bc_palette["green"])
)

# ---- Shared theme: theme_classic + top legend + bold facet strips.
theme_cae <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "top",
      strip.text.x    = ggplot2::element_text(face = "bold"),
      strip.text.y    = ggplot2::element_text(face = "bold"))
}

# ---- pretty()-based axis ticks (define once, reuse).
number_ticks <- function(n) function(limits) pretty(limits, n)

# ---- Probe whether cairo_pdf actually works at runtime. capabilities("cairo")
# can report TRUE while the cairo DLL fails to load (missing X11 libs on some
# macOS installs), so we test it for real once and cache the result.
.cae_cairo_ok <- local({
  cached <- NA
  function() {
    if (!is.na(cached)) return(cached)
    ok <- FALSE
    if (isTRUE(capabilities("cairo"))) {
      tf <- tempfile(fileext = ".pdf")
      ok <- tryCatch(
        withCallingHandlers({
          grDevices::cairo_pdf(tf)
          try(grDevices::dev.off(), silent = TRUE)
          isTRUE(file.exists(tf) && file.info(tf)$size > 0)
        },
        warning = function(w) invokeRestart("muffleWarning")),
        error = function(e) FALSE)
      if (file.exists(tf)) unlink(tf)
    }
    cached <<- isTRUE(ok)
    if (!cached)
      message("[plot_theme] cairo_pdf unavailable; using base 'pdf' device for vector output")
    cached
  }
})

# ---- Locate a Ghostscript binary for font embedding, once.
.cae_gs <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    cand <- c(Sys.getenv("R_GSCMD", ""), Sys.which("gs"),
              "/opt/homebrew/bin/gs", "/usr/local/bin/gs", "/usr/bin/gs")
    cand <- cand[nzchar(cand)]
    hit <- cand[file.exists(cand)]
    cached <<- if (length(hit)) hit[1] else NA_character_
    cached
  }
})

# ---- Embed all fonts into a vector PDF, in place.
# The base 'pdf' device references the standard-14 fonts (Helvetica) by name
# without embedding them. The journal's artwork guidelines require embedded
# fonts, and `pdffonts` shows `emb: no` for the base device. Ghostscript rewrites the file with subsetted Type 1C
# outlines and leaves the vector geometry editable. cairo_pdf embeds fonts on its
# own, so this is a no-op path when cairo is available.
.embed_pdf_fonts <- function(path) {
  gs <- .cae_gs()
  if (is.na(gs)) {
    warning("[plot_theme] Ghostscript not found; fonts NOT embedded in ", basename(path),
            ". Install ghostscript or set R_GSCMD.", call. = FALSE)
    return(invisible(FALSE))
  }
  tmp <- paste0(path, ".embed.tmp")
  status <- system2(gs,
    c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
      "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dPDFSETTINGS=/prepress",
      "-dCompatibilityLevel=1.4",
      paste0("-sOutputFile=", shQuote(tmp)), shQuote(path)),
    stdout = FALSE, stderr = FALSE)
  if (status == 0 && file.exists(tmp) && file.info(tmp)$size > 0) {
    file.rename(tmp, path)
    invisible(TRUE)
  } else {
    if (file.exists(tmp)) unlink(tmp)
    warning("[plot_theme] Ghostscript failed to embed fonts in ", basename(path),
            call. = FALSE)
    invisible(FALSE)
  }
}

# ---- Save a figure as BOTH 600-dpi PNG and an editable vector PDF.
# Lancet requires editable vector artwork with embedded fonts; every figure is
# written twice, and the PDF twin is post-processed to embed its fonts.
save_fig <- function(plot, name, width, height, output_dir = CFG$dir_fig) {
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(output_dir, paste0(name, ".png"))
  pdf_path <- file.path(output_dir, paste0(name, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 600,
                  bg = "white")
  cairo <- .cae_cairo_ok()
  pdf_device <- if (cairo) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  device = pdf_device, bg = "white")
  if (!cairo) .embed_pdf_fonts(pdf_path)
  invisible(c(png_path, pdf_path))
}
