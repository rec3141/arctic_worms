#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Data-prep (run once, needs the raw DADA2 tables which are NOT in this repo):
# aggregate the per-sample Platyhelminthes ASV table into a per-SITE summary of
# mean relative abundance, so only site-level data is distributed.
# Produces  data/asv_site_summary.tsv  and  data/locations.tsv .
# ---------------------------------------------------------------------------
RAW <- "/Users/ericcollins/Downloads/project_Acoela/work"                # raw inputs (not committed)
OUT <- "/Users/ericcollins/Downloads/project_Acoela/arctic_worms/data"

## locations (GPS) + canonical site names -----------------------------------
rl <- readLines(file.path(RAW, "Readme.txt")); rl <- rl[grep("^locations", rl)[1]:length(rl)][-1]
rl <- rl[nzchar(trimws(rl))]; num <- "-?[0-9]+[.][0-9]+"
locs <- do.call(rbind, lapply(rl, function(x){ x <- gsub("−","-",x)
  m <- regmatches(x, regexpr(paste0(num,",\\s*",num), x)); if(!length(m)) return(NULL)
  ll <- as.numeric(strsplit(m, ",")[[1]]); data.frame(name=trimws(sub(m,"",x,fixed=TRUE)), lat=ll[1], lon=ll[2]) }))
locs <- locs[!locs$name %in% c("CSTC","SRR","TRAVERSE"), ]     # off-map / skipped
S2 <- locs$name[grepl("S2-I-R2",locs$name)]; MO <- locs$name[grepl("MO-D1",locs$name)]; CH <- locs$name[grepl("CH5A",locs$name)]
locs$region <- locs$name
locs$region[locs$name==S2] <- "UTQ"; locs$region[locs$name==MO] <- "MinMOS"
locs$region[locs$name==CH] <- "CORALHARBOUR"; locs$region[locs$name=="Alex"] <- "SITKA"
region_of <- function(sid){ s <- toupper(sid)
  if (grepl("^S[0-9]+-[IK]-",s)) return("UTQ"); if (grepl("MO-D|^EL-MO|^C2",s)) return("MinMOS")
  if (grepl("-CH[0-9]+",s)) return("CORALHARBOUR"); if (grepl("^ALEX",s)) return("SITKA")
  for (nm in locs$region[!grepl("^like",locs$name)]) if (grepl(paste0("^",toupper(nm)),s)) return(nm)
  NA_character_ }
write.table(locs[,c("region","lat","lon")], file.path(OUT,"locations.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

## per-site mean relative abundance per ASV ---------------------------------
seqids <- read.table(file.path(RAW,"analysis_acoela/seqids.txt"), sep="\t", header=TRUE)
nt <- as.matrix(read.table(file.path(RAW,"analysis_acoela/normtab_platy2.tsv"), sep="\t",
                           header=TRUE, row.names=1, check.names=FALSE, quote="\""))
fa <- readLines(file.path(RAW,"analysis_acoela/Platyhelminthes2.fasta")); h <- grep("^>",fa)
ids <- sub("^>(ASV[0-9]+)_.*","\\1",fa[h])
seqs <- vapply(seq_along(h), function(i) paste(fa[(h[i]+1):(if(i<length(h)) h[i+1]-1 else length(fa))],collapse=""), character(1))
colnames(nt) <- setNames(ids,seqs)[colnames(nt)]
reg <- vapply(seqids$sample_id[match(rownames(nt), seqids$seqid)], region_of, character(1))
keep <- !is.na(reg); nt <- nt[keep,]; reg <- reg[keep]

sites <- sort(unique(reg))
M <- t(sapply(sites, function(r) colMeans(nt[reg==r, , drop=FALSE])))     # region x ASV mean RA
n_samples <- as.integer(table(reg)[sites])
out <- data.frame(region=sites, n_samples=n_samples, round(as.data.frame(M), 6), check.names=FALSE)
write.table(out, file.path(OUT,"asv_site_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat("wrote asv_site_summary.tsv:", nrow(out), "sites x", ncol(M), "ASVs;  locations.tsv:", nrow(locs), "sites\n")
