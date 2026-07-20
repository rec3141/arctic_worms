#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Full publication placement tree for Arctic flatworm 18S phylogenetic
# placement.  Reference tree (grey, unlabelled) provides context; the KBD
# cultured isolates, placed ASV amplicons, and the outgroup are highlighted.
# Rooted on the outgroup; ASVs filtered by placement pendant length (<0.2).
#
# Usage:  Rscript plot_placement_tree.R <group>
#   group in {acoela, rhabdocoela, macrostomorpha}  (default: acoela)
# ---------------------------------------------------------------------------
suppressMessages({library(ape)})

args  <- commandArgs(trailingOnly = TRUE)
group <- ifelse(length(args) >= 1, args[1], "acoela")

cfg <- list(
  acoela         = list(dir = "../data/trees/acoela",         title = "Acoela"),
  rhabdocoela    = list(dir = "../data/trees/rhabdocoela",    title = "Rhabdocoela"),
  macrostomorpha = list(dir = "../data/trees/macrostomorpha", title = "Macrostomorpha")
)[[group]]
stopifnot(!is.null(cfg))

PENDANT_MAX <- 0.2

kbd_species <- c(
  "15"="Baltalimania naatanaŋŋuakulluraq","16"="Baltalimania naatanaŋŋuakulluraq",
  "24"="Pseudohaplogonaria kakilisaliuralik","21"="red acoel (undescribed)",
  "11"="Baicalellia sp.","19"="Promesostoma sp.","04"="undescribed sp.",
  "07"="sea-ice acoel","09"="sediment acoel","13"="sediment acoel","06"="undescribed sp.")

# --- read, pendant-filter, root on outgroup --------------------------------
tr <- read.tree(file.path(cfg$dir, "tog.nwk"))
p  <- read.csv(file.path(cfg$dir, "placements.csv"), stringsAsFactors = FALSE)
best <- do.call(rbind, lapply(split(p, p$name), function(d) d[which.max(d$like_weight_ratio), ]))
drop_asv <- intersect(best$name[best$pendant_length >= PENDANT_MAX], tr$tip.label)
if (length(drop_asv)) tr <- drop.tip(tr, drop_asv)

# drop uchime-flagged chimeras (pendant length misses divergent-parent chimeras)
if (file.exists("../data/chimera_results.txt")) {
  cr <- read.table("../data/chimera_results.txt", sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
  chim <- sub("_.*", "", cr$V2[cr[[ncol(cr)]] == "Y"])
  drop_chim <- tr$tip.label[sub("_.*", "", tr$tip.label) %in% chim]
  if (length(drop_chim)) tr <- drop.tip(tr, drop_chim)
}

# root on the edge between the outgroup and the ingroup (ingroup-MRCA branch)
og_tips <- grep("OUTGRP", tr$tip.label, value = TRUE)
ingroup <- setdiff(tr$tip.label, og_tips)
if (length(og_tips) >= 1 && length(ingroup) >= 2)
  tr <- tryCatch(root(tr, node = getMRCA(tr, ingroup), resolve.root = TRUE),
                 error = function(e) root(tr, outgroup = og_tips[1], resolve.root = TRUE))
tr <- ladderize(tr)
tips <- tr$tip.label

is_og  <- grepl("OUTGRP", tips)
is_kbd <- grepl("KBD", tips) & !is_og
is_asv <- grepl("ASV", tips) & !is_og
is_ref <- !is_kbd & !is_asv & !is_og

# --- pretty labels ---------------------------------------------------------
lab <- character(length(tips))
kbd_num <- sub(".*KBD_?([0-9]+).*", "\\1", tips); kbd_sp <- kbd_species[kbd_num]
lab[is_kbd] <- ifelse(is.na(kbd_sp[is_kbd]), paste0("KBD-", kbd_num[is_kbd]),
                      paste0("KBD-", kbd_num[is_kbd], "  ", kbd_sp[is_kbd]))
lab[is_asv] <- sub("^(ASV[0-9]+)_.*", "\\1", tips[is_asv])   # ASV number only
og_lab <- sub("_OUTGRP.*$", "", sub("^_R_", "", tips))
og_lab <- sub("_([A-Z]{2}[0-9]{4,})$", " \\1", og_lab)
lab[is_og] <- paste0(gsub("_", " ", og_lab[is_og]), "  (outgroup)")
lab[is_ref] <- ""
tr$tip.label <- lab

# --- colours & sizes -------------------------------------------------------
col_ref <- "grey78"; col_kbd <- "#C1272D"; col_asv <- "#1F6FB2"; col_og <- "#2E7D5B"
tipcol <- ifelse(is_kbd, col_kbd, ifelse(is_asv, col_asv, ifelse(is_og, col_og, col_ref)))
tipcex <- ifelse(is_kbd, 0.68, ifelse(is_asv, 0.55, ifelse(is_og, 0.68, 0.0001)))
tipfont<- ifelse(is_kbd, 2, ifelse(is_asv, 3, ifelse(is_og, 3, 1)))
ptcol  <- tipcol
ptcex  <- ifelse(is_kbd, 1.5, ifelse(is_asv, 0.9, ifelse(is_og, 1.3, 0.35)))

# --- draw ------------------------------------------------------------------
n <- length(tips); h <- max(8, n * 0.115); depth <- max(node.depth.edgelength(tr))
cairo_pdf(file.path("../figures", paste0(group, "_placement_tree.pdf")), width = 12, height = h)
par(mar = c(1, 1, 3, 1), xpd = NA)
plot.phylo(tr, type = "phylogram", cex = tipcex, font = tipfont, tip.color = tipcol,
           label.offset = depth * 0.006, edge.color = "grey45", edge.width = 0.6,
           x.lim = c(0, depth * 1.45), y.lim = c(-1.5, n + 0.5))
tiplabels(pch = 19, col = ptcol, cex = ptcex)
add.scale.bar(x = depth * 0.02, y = -0.7, length = 0.05, cex = 0.8, lwd = 1.4, col = "grey20")
title(main = paste0(cfg$title, " — 18S phylogenetic placement"), cex.main = 1.4, font.main = 1)
legend("topright", inset = c(0.02, 0.02), bty = "n", pt.cex = c(1.5, 0.9, 1.3, 0.6),
       pch = 19, col = c(col_kbd, col_asv, col_og, col_ref), text.font = c(2, 3, 3, 1), cex = 0.9,
       legend = c("KBD cultured isolate", "ASV amplicon (placed)", "Outgroup", "Reference 18S"))
invisible(dev.off())
cat("wrote out/", group, "_placement_tree.pdf  (", n, " tips )\n", sep = "")
