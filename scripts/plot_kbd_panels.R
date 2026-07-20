#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Composable placement-figure panels (KBD overview + zoom panels) for a group.
#   Usage:  Rscript plot_kbd_panels.R <group>   (default: acoela)
# Outputs one PDF per panel in ./out so they can be laid out together.
# ---------------------------------------------------------------------------
suppressMessages(library(ape))

grp_arg <- commandArgs(trailingOnly = TRUE)
group   <- if (length(grp_arg) >= 1) grp_arg[1] else "acoela"

CFG <- list(
  acoela = list(dir = "../data/trees/acoela",
    title = "Acoela — placement of KBD isolates",
    panels = list(
      list(file="acoela_panel_g1_KBD24_KBD13.pdf", nums=c("24","13"), up=0,
           main="KBD-24 (Pseudohaplogonaria kakilisaliuralik) + KBD-13"),
      list(file="acoela_panel_g2_KBD09.pdf",       nums=c("09"),      up=3,
           main="KBD-09 (Haploposthia neighbourhood)"),
      list(file="acoela_panel_g3_KBD15_KBD16.pdf", nums=c("15","16"), up=2,
           main="KBD-15 + KBD-16 (Baltalimania naatanaŋŋuakulluraq)"),
      list(file="acoela_panel_g4_KBD07_KBD21.pdf", nums=c("07","21"), up=3,
           main="KBD-07 + KBD-21 (Haplogonaria radiation)"))),
  rhabdocoela = list(dir = "../data/trees/rhabdocoela",
    title = "Rhabdocoela — placement of KBD isolates",
    panels = list(
      list(file="rhabdocoela_panel_r1_Baicalellia.pdf", nums=c("04","06","11"), up=4,
           main="KBD-04 + KBD-06 + KBD-11 (Baicalellia)"),
      list(file="rhabdocoela_panel_r2_KBD19.pdf",       nums=c("19"),           up=6,
           main="KBD-19 (Promesostoma)")))
)
cfg <- CFG[[group]]; stopifnot(!is.null(cfg)); outdir <- cfg$dir

tr  <- read.tree(file.path(outdir, "tog.nwk"))   # all placements (best position each)

# --- keep only confidently-anchored ASVs: pendant length < 0.2 -------------
# (removes amplicons placed on long pendant branches, but retains ASVs even
#  when placement mass splits among near-identical tips)
PENDANT_MAX <- 0.2
p <- read.csv(file.path(outdir, "placements.csv"), stringsAsFactors = FALSE)
best <- do.call(rbind, lapply(split(p, p$name),
          function(d) d[which.max(d$like_weight_ratio), ]))
drop_asv <- intersect(best$name[best$pendant_length >= PENDANT_MAX], tr$tip.label)
if (length(drop_asv)) tr <- drop.tip(tr, drop_asv)

# --- drop uchime-flagged chimeras (belt-and-suspenders on DADA2) ------------
# pendant length does NOT catch chimeras between divergent parents (e.g. ASV149
# = Baicalellia x Promesostoma, pendant 0.13), so remove them explicitly.
if (file.exists("../data/chimera_results.txt")) {
  cr <- read.table("../data/chimera_results.txt", sep = "\t", fill = TRUE,
                    stringsAsFactors = FALSE, quote = "")
  chim <- sub("_.*", "", cr$V2[cr[[ncol(cr)]] == "Y"])          # ASV numbers
  drop_chim <- tr$tip.label[sub("_.*", "", tr$tip.label) %in% chim]
  if (length(drop_chim)) tr <- drop.tip(tr, drop_chim)
}

# Root on the edge between the outgroup and the ingroup: place the root on the
# branch subtending the ingroup MRCA, so the (possibly paraphyletic) outgroup
# grade sits on one side and the whole ingroup on the other.
og_tips  <- grep("OUTGRP", tr$tip.label, value = TRUE)
ingroup  <- setdiff(tr$tip.label, og_tips)
if (length(og_tips) >= 1 && length(ingroup) >= 2) {
  tr <- tryCatch(root(tr, node = getMRCA(tr, ingroup), resolve.root = TRUE),
                 error = function(e) root(tr, outgroup = og_tips[1], resolve.root = TRUE))
}
tr  <- ladderize(tr)
tips <- tr$tip.label

col_ref <- "grey78"; col_kbd <- "#C1272D"; col_asv <- "#1F6FB2"
col_og <- "#2E7D5B"; col_reftxt <- "grey25"

kbd_species <- c(
  "15" = "Baltalimania naatanaŋŋuakulluraq",
  "16" = "Baltalimania naatanaŋŋuakulluraq",
  "24" = "Pseudohaplogonaria kakilisaliuralik",
  "21" = "red acoel (undescribed)",
  "07" = "sea-ice acoel", "09" = "sediment acoel", "13" = "sediment acoel",
  "11" = "Baicalellia sp.", "19" = "Promesostoma sp.",
  "04" = "undescribed sp.", "06" = "undescribed sp."
)

cat_of  <- function(x) ifelse(grepl("OUTGRP",x),"OG",
                       ifelse(grepl("KBD",x),"KBD",
                       ifelse(grepl("ASV",x),"ASV","ref")))
kbd_num <- function(x) sub(".*KBD_?([0-9]+).*", "\\1", x)
fmt_og  <- function(x){                       # Nemertoderma_westbladi_AF327726_OUTGRP_18S
  x <- sub("_OUTGRP.*$", "", x); x <- sub("^_R_", "", x)
  x <- sub("_([A-Z]{2}[0-9]{4,})$", " \\1", x); gsub("_", " ", x)
}

# --- pretty labels ---------------------------------------------------------
fmt_ref <- function(x){                       # Acoela_Genus_sp_ID_xxx_18S -> "Genus sp xxx"
  x <- sub("_18S$", "", x); x <- sub("^Acoela_", "", x)
  x <- gsub("_(ID|GB)_", " ", x); gsub("_", " ", x)
}
fmt_asv <- function(x) sub("^(ASV[0-9]+)_.*", "\\1", x)   # ASV number only
# (the "_Genus_species" in ASV names is an automated DADA2/PR2 assignment and
#  is unreliable at species level, so we show only the ASV id)
fmt_kbd <- function(x){
  n <- kbd_num(x); sp <- kbd_species[n]
  ifelse(is.na(sp), paste0("KBD-", n), paste0("KBD-", n, "  ", sp))
}

# node = MRCA of the KBD tips, then walk up `up` parent edges
clade_node <- function(labs, up){
  idx <- match(labs, tips)
  anc <- if(length(idx) >= 2) getMRCA(tr, labs) else tr$edge[tr$edge[,2]==idx, 1]
  for(i in seq_len(up)){
    par <- tr$edge[tr$edge[,2]==anc, 1]; if(length(par)==0) break; anc <- par
  }
  anc
}

draw_tree <- function(t, file, main, show_ref_labels, height, width=9){
  tl  <- t$tip.label; ct <- cat_of(tl)
  lab <- ifelse(ct=="KBD", fmt_kbd(tl),
         ifelse(ct=="ASV", fmt_asv(tl),
         ifelse(ct=="OG",  paste0(fmt_og(tl), "  (outgroup)"),
                if(show_ref_labels) fmt_ref(tl) else "")))
  # extra leading space on KBD labels so the red dot doesn't touch the text
  lab <- ifelse(ct=="KBD", paste0("   ", lab), lab)
  t$tip.label <- lab
  tcol <- ifelse(ct=="KBD", col_kbd, ifelse(ct=="ASV", col_asv, ifelse(ct=="OG", col_og, col_reftxt)))
  tcex <- ifelse(ct=="KBD", 0.85, ifelse(ct=="ASV", 0.72, ifelse(ct=="OG", 0.8, 0.66)))
  tfnt <- ifelse(ct=="KBD", 2,    ifelse(ct=="ASV", 1,    3))
  pcol <- ifelse(ct=="KBD", col_kbd, ifelse(ct=="ASV", col_asv, ifelse(ct=="OG", col_og, col_ref)))
  pcex <- ifelse(ct=="KBD", 1.7, ifelse(ct=="ASV", 1.0, ifelse(ct=="OG", 1.3, 0.7)))
  depth <- max(node.depth.edgelength(t)); n <- length(tl); off <- depth * 0.012
  # size the x extent so every tip label fits: use real tip x-positions + label widths (inches)
  tf <- tempfile(fileext=".pdf"); pdf(tf, width=width, height=height); par(mar=c(1,1,3,1))
  plot.phylo(t, cex=tcex, font=tfnt, root.edge=TRUE, label.offset=off, plot=FALSE, x.lim=c(0, depth*2))
  xx <- get("last_plot.phylo", envir=.PlotPhyloEnv)$xx[seq_len(n)]
  pw <- par("pin")[1]
  lw <- vapply(seq_len(n), function(i) strwidth(t$tip.label[i], units="inches", cex=tcex[i], font=tfnt[i]), numeric(1))
  dev.off(); unlink(tf)
  xr <- 1.02 * max((xx + off) * pw / pmax(pw - lw - 0.05, 0.15))   # per-tip requirement
  cairo_pdf(file.path("../figures", file), width=width, height=height)
  par(mar=c(1,1,3,1), xpd=NA)
  plot.phylo(t, cex=tcex, font=tfnt, tip.color=tcol, root.edge=TRUE,
             label.offset=off, edge.color="grey45", edge.width=0.8,
             x.lim=c(0, xr), y.lim=c(-1.5, n + 0.5))
  tiplabels(pch=19, col=pcol, cex=pcex)
  add.scale.bar(x=depth*0.02, y=-0.7, length=0.05, cex=0.75, lwd=1.3, col="grey20")
  title(main=main, cex.main=1.25, font.main=1)
  invisible(dev.off())
  cat("wrote out/", file, "  (", length(tl), " tips)\n", sep="")
}

# --- KBD overview, three versions for comparison ---------------------------
#   asis      : grey context tree, only KBD + outgroup labelled (leader lines)
#   reflabels : as-is plus every reference tip labelled (grey)
#   pruned    : tree pruned to KBD + outgroup leaves only
draw_overview <- function(mode) {
  t <- tr
  if (mode == "pruned") {
    t <- keep.tip(t, tips[cat_of(tips) %in% c("KBD", "OG")])
  } else {                                # trim ASV hits from the context tree
    asv <- tips[cat_of(tips) == "ASV"]
    if (length(asv)) t <- drop.tip(t, asv)
  }
  tl <- t$tip.label; ct <- cat_of(tl)
  depth <- max(node.depth.edgelength(t)); n <- length(tl)
  H <- if (mode == "pruned") max(3.5, n * 0.32) else 11
  suffix <- c(asis="", reflabels="_reflabels", pruned="_pruned")[mode]
  cairo_pdf(sprintf("../figures/%s_KBD_overview%s.pdf", group, suffix), width=8.5, height=H)
  par(mar=c(1,1,3,1), xpd=NA)

  if (mode == "pruned") {
    lab <- ifelse(ct=="KBD", paste0("  ", fmt_kbd(tl)),
                  paste0("  ", fmt_og(tl), " (outgroup)"))
    t$tip.label <- lab
    tcol <- ifelse(ct=="KBD", col_kbd, col_og)
    plot.phylo(t, cex=0.9, font=ifelse(ct=="OG",3,2), tip.color=tcol, root.edge=TRUE,
               label.offset=depth*0.02, edge.color="grey45", edge.width=1.1,
               x.lim=c(0, depth*1.7), y.lim=c(-1.2, n+0.5))
    tiplabels(pch=19, col=tcol, cex=ifelse(ct=="KBD",1.9,1.3))
    add.scale.bar(x=depth*0.02, y=-0.6, length=0.05, cex=0.75, lwd=1.3, col="grey20")
    title(main=cfg$title, cex.main=1.3, font.main=1)
    invisible(dev.off()); cat("wrote out/", group, "_KBD_overview", suffix, ".pdf  (", n, " tips)\n", sep=""); return(invisible())
  }

  plot.phylo(t, show.tip.label=FALSE, edge.color="grey55", edge.width=0.6,
             x.lim=c(0, depth*1.9), y.lim=c(-2, n + 0.5))
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv); xx <- pp$xx[1:n]; yy <- pp$yy[1:n]
  if (mode == "reflabels") {                     # grey label on every reference tip
    ri <- which(ct=="ref")
    # shrink font until adjacent labels no longer overlap (tips are 1 user-unit apart)
    cex_ref <- min(0.42, 0.92 / strheight("Ag", cex = 1, units = "user"))
    text(xx[ri]+depth*0.006, yy[ri], fmt_ref(tl[ri]), col=col_reftxt, font=3, cex=cex_ref, adj=0)
  }
  points(xx[ct=="KBD"], yy[ct=="KBD"], pch=19, col=col_kbd, cex=1.9)
  points(xx[ct=="OG"],  yy[ct=="OG"],  pch=19, col=col_og,  cex=1.3)
  li  <- which(ct %in% c("KBD","OG")); li <- li[order(yy[li])]
  ly  <- yy[li]; ming <- n*0.02
  for (k in seq_along(ly)[-1]) if (ly[k]-ly[k-1] < ming) ly[k] <- ly[k-1]+ming
  for (k in seq_along(li)) {
    i <- li[k]
    txt <- if (ct[i]=="OG") paste0(fmt_og(tl[i])," (outgroup)") else fmt_kbd(tl[i])
    col <- if (ct[i]=="OG") col_og else col_kbd
    segments(xx[i]+depth*0.01, yy[i], xx[i]+depth*0.03, ly[k], col=col, lwd=0.5)
    text(xx[i]+depth*0.035, ly[k], txt, col=col, font=if(ct[i]=="OG")3 else 2, cex=0.8, adj=0)
  }
  add.scale.bar(x=depth*0.02, y=-1, length=0.05, cex=0.75, lwd=1.3, col="grey20")
  title(main=cfg$title, cex.main=1.3, font.main=1)
  invisible(dev.off())
  cat("wrote out/", group, "_KBD_overview", suffix, ".pdf  (", n, " tips)\n", sep="")
}
for (m in c("asis","reflabels","pruned")) draw_overview(m)

# --- zoom panels (from group config) ---------------------------------------
for(p in cfg$panels){
  labs <- unlist(lapply(p$nums, function(n) grep(paste0("KBD_?", n, "_"), tips, value=TRUE)))
  cl   <- extract.clade(tr, clade_node(labs, p$up), root.edge = 1)  # keep stem to root
  h    <- max(3.2, length(cl$tip.label) * 0.22)
  draw_tree(cl, p$file, p$main, show_ref_labels=TRUE, height=h)
}
