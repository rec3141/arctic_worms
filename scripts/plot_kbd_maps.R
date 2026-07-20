#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# KBD-isolate Arctic distribution maps (run from scripts/):
#   <isolate>_map_presence.pdf  - per-isolate presence / absence by site
#   <group>_map_pie.pdf         - per-group pie: cumulative rel. abundance / isolate
#   <group>_map_richness.pdf    - per-group pie: # distinct ASVs / isolate
# Uses the SITE-LEVEL summary (../data/asv_site_summary.tsv), not raw samples.
# KBD->ASV assignment from the placement trees (nearest isolate, <2% divergence).
# ---------------------------------------------------------------------------
suppressMessages({library(ape); library(sf); library(ggplot2); library(scatterpie)
                  library(rnaturalearth)})
FIG <- "../figures"; ASV_MAX_DIV <- 0.02

## --- 1) sites + per-site ASV mean relative abundance ----------------------
locs <- read.table("../data/locations.tsv", sep = "\t", header = TRUE, check.names = FALSE)
locs$name <- locs$region; locs$label <- locs$region
S  <- read.table("../data/asv_site_summary.tsv", sep = "\t", header = TRUE, check.names = FALSE)
M0 <- as.matrix(S[, setdiff(colnames(S), c("region","n_samples"))]); rownames(M0) <- S$region
iso_ab   <- function(ids) { ids <- intersect(ids, colnames(M0)); if(!length(ids)) return(NULL)
                            rowSums(M0[, ids, drop = FALSE]) }                 # region -> summed mean RA
iso_rich <- function(ids) { ids <- intersect(ids, colnames(M0)); if(!length(ids)) return(setNames(rep(0,nrow(M0)),rownames(M0)))
                            rowSums(M0[, ids, drop = FALSE] > 0) }             # region -> # ASVs detected

## --- 2) KBD -> ASV per group (from placement trees) -----------------------
kbd_asv <- list()
for (g in c("acoela","rhabdocoela","macrostomorpha")) {
  f <- paste0("../data/trees/", g, "/tog.nwk"); if (!file.exists(f)) next
  t <- read.tree(f); tips <- t$tip.label; kbd <- grep("KBD", tips, value = TRUE); if (!length(kbd)) next
  D <- cophenetic(t)
  for (a in grep("ASV", tips, value = TRUE)) {
    dk <- D[a, kbd]; if (min(dk) > ASV_MAX_DIV) next
    num <- sub(".*KBD_?([0-9]+).*", "\\1", kbd[which.min(dk)])
    kbd_asv[[g]][[num]] <- union(kbd_asv[[g]][[num]], sub("^(ASV[0-9]+)_.*", "\\1", sub("^_R_", "", a)))
  }
  if (g == "acoela" && !is.null(kbd_asv[[g]][["16"]])) {
    kbd_asv[[g]][["15"]] <- union(kbd_asv[[g]][["15"]], kbd_asv[[g]][["16"]]); kbd_asv[[g]][["16"]] <- NULL
  }
}
disp <- c("15"="Baltalimania naatanaŋŋuakulluraq","07"="KBD-07 (Haplogonaria)",
          "21"="KBD-21 (red acoel)","09"="KBD-09 (Haploposthia)","19"="KBD-19 (Promesostoma)")
pal  <- c("15"="#C1272D","07"="#1F6FB2","09"="#2E7D5B","21"="#E39A2C","19"="#7A4FA3")

## --- 3) base map ----------------------------------------------------------
crs_polar <- "+proj=laea +lat_0=90 +lon_0=-150 +datum=WGS84 +units=m"
suppressMessages(sf_use_s2(FALSE))
land <- st_transform(st_crop(ne_countries(scale="medium", returnclass="sf"),
                             xmin=-180, xmax=180, ymin=40, ymax=90), crs_polar)
xy <- st_coordinates(st_transform(st_as_sf(locs, coords = c("lon","lat"), crs = 4326), crs_polar))
locs$X <- xy[,1]; locs$Y <- xy[,2]
lim <- list(x = range(locs$X) + c(-4e5, 9e5), y = range(locs$Y) + c(-4e5, 4e5))
base_map <- function()
  ggplot() +
    geom_sf(data = land, fill = "grey92", colour = "grey70", linewidth = 0.2) +
    coord_sf(crs = crs_polar, xlim = lim$x, ylim = lim$y, expand = FALSE) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 11) +
    theme(panel.background = element_rect(fill = "#eef3f7"),
          panel.grid = element_line(colour = "grey88"),
          axis.text = element_blank(), axis.ticks = element_blank())

## --- 4a) per-isolate presence / absence maps ------------------------------
for (g in names(kbd_asv)) for (num in names(kbd_asv[[g]])) {
  ab <- iso_ab(kbd_asv[[g]][[num]]); if (is.null(ab)) next
  st <- merge(locs, data.frame(name = names(ab), present = as.integer(ab > 0)), by = "name")
  sp <- ifelse(is.na(disp[num]), paste0("KBD-", num), disp[num]); safe <- gsub("[^A-Za-z0-9]+","_", sp)
  p <- base_map() +
    geom_point(data = st, aes(X, Y, fill = factor(present)), shape = 21, size = 4, stroke = 0.5, colour = "grey20") +
    ggrepel::geom_text_repel(data = st, aes(X, Y, label = label), size = 2.6, colour = "grey25", max.overlaps = 20) +
    scale_fill_manual(values = c("0"="white","1"=unname(pal[num])), labels = c("absent","present"), name = NULL) +
    labs(title = paste0(sp, " — presence / absence"))
  ggsave(file.path(FIG, paste0(safe, "_map_presence.pdf")), p, width = 7.5, height = 6.5, device = cairo_pdf)
}

## --- 4b) per-group pie maps (cumulative RA, and ASV richness) --------------
draw_pie_map <- function(g, nums, M, title, suffix) {
  M[is.na(M)] <- 0; colnames(M) <- nums
  w <- cbind(locs, as.data.frame(M)); w$total <- rowSums(M); w <- w[w$total > 0, ]
  w$r <- 1.1e5 + 2.2e5 * sqrt(w$total / max(w$total))
  ROT <- pi/2; MIN_SLICE <- 0.03
  wd <- do.call(rbind, lapply(order(-w$total), function(ri) {
    row <- w[ri, ]; vals <- as.numeric(row[nums]); vals[vals / sum(vals) < MIN_SLICE] <- 0
    tot <- sum(vals); nz <- which(vals > 0)
    if (length(nz) == 1) { th <- seq(0, 2*pi, length.out = 100)
      return(data.frame(X = row$X + row$r*sin(th), Y = row$Y + row$r*cos(th),
                        iso = nums[nz], sid = paste0(row$name, "_", nums[nz]), ord = ri)) }
    a <- ROT; segs <- list()
    for (j in nz) { a0 <- a; a1 <- a + vals[j]/tot * 2*pi; a <- a1
      th <- seq(a0, a1, length.out = max(2, ceiling(vals[j]/tot * 80)))
      segs[[length(segs)+1]] <- data.frame(X = c(row$X, row$X + row$r*sin(th)), Y = c(row$Y, row$Y + row$r*cos(th)),
        iso = nums[j], sid = paste0(row$name, "_", nums[j]), ord = ri) }
    do.call(rbind, segs)
  }))
  wd <- wd[order(-w$total[match(sub("_.*","", wd$sid), w$name)]), ]; wd$sid <- factor(wd$sid, levels = unique(wd$sid))
  none <- locs[!locs$name %in% w$name, ]
  p <- base_map() +
    geom_point(data = none, aes(X, Y), shape = 4, size = 2.2, stroke = 0.9, colour = "grey55") +
    geom_polygon(data = wd, aes(X, Y, group = sid, fill = iso), colour = "grey30", linewidth = 0.2) +
    ggrepel::geom_text_repel(data = locs, aes(X, Y, label = label), size = 2.6, colour = "grey25",
                             max.overlaps = 20, point.padding = 8) +
    scale_fill_manual(values = pal[nums], labels = ifelse(is.na(disp[nums]), paste0("KBD-",nums), disp[nums]),
                      name = "KBD isolate", breaks = nums) +
    labs(title = title)
  ggsave(file.path(FIG, paste0(g, "_map_", suffix, ".pdf")), p, width = 8.5, height = 6.5, device = cairo_pdf)
  cat("wrote", g, suffix, "map (", length(nums), "isolates, ", nrow(w), "sites )\n")
}
for (g in names(kbd_asv)) {
  nums <- names(kbd_asv[[g]]); if (!length(nums)) next
  Mab <- sapply(nums, function(n) iso_ab(kbd_asv[[g]][[n]])[locs$name])
  Mri <- sapply(nums, function(n) iso_rich(kbd_asv[[g]][[n]])[locs$name])
  draw_pie_map(g, nums, Mab, paste0(tools::toTitleCase(g), " — cumulative relative abundance"), "pie")
  draw_pie_map(g, nums, Mri, paste0(tools::toTitleCase(g), " — ASV richness"), "richness")
}
