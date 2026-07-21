# Render figure text as vector outlines (via showtext) so PDFs display
# identically in every viewer — in particular macOS Preview, which mangles the
# subsetted synthetic-oblique fonts that cairo_pdf embeds by default.
# Falls back to the device default font if showtext / a suitable font is absent.
FIGFAM <- ""   # "" = device default
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  cands <- list(
    c("/System/Library/Fonts/Supplemental/Arial.ttf",
      "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
      "/System/Library/Fonts/Supplemental/Arial Italic.ttf"),
    c("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
      "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
      "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf"),
    c("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
      "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
      "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf"))
  for (f in cands) if (all(file.exists(f))) {
    sysfonts::font_add("fig", regular = f[1], bold = f[2], italic = f[3])
    showtext::showtext_auto(); showtext::showtext_opts(dpi = 72); FIGFAM <- "fig"; break
  }
}
