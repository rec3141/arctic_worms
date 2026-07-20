#!/bin/bash
# ---------------------------------------------------------------------------
# Build a rooted 18S phylogenetic-placement tree for one taxon group:
# dereplicate + align references(+isolates)+outgroup, infer an ML tree, add
# SH-like support, then place the amplicon ASVs onto it and filter.
#
#   run_pipeline.sh <outdir> <references.fasta> <outgroups.fasta> \
#                   <amplicons.fasta> [isolates.fasta] [nsearch]
#
# Requires on PATH: vsearch mafft raxmlHPC-PTHREADS-AVX taxit pplacer guppy gappa
# (taxtastic / pplacer / gappa are typically each in their own env; activate
#  them however your system provides them before running.)
# ---------------------------------------------------------------------------
set -euo pipefail
OUTDIR=$1; REFS=$2; OG=$3; SHORT=$4; ISO=${5:-/dev/null}; NSEARCH=${6:-1000}

# pick whichever threaded RAxML binary is installed (bioconda ships -AVX2)
RAXML=$(command -v raxmlHPC-PTHREADS-AVX2 || command -v raxmlHPC-PTHREADS-AVX \
     || command -v raxmlHPC-PTHREADS-SSE3 || command -v raxmlHPC-PTHREADS \
     || command -v raxmlHPC) || { echo "no RAxML binary found" >&2; exit 1; }

# resolve inputs to absolute paths (we cd into OUTDIR below)
REFS=$(realpath "$REFS"); OG=$(realpath "$OG"); SHORT=$(realpath "$SHORT"); ISO=$(realpath "$ISO")

mkdir -p "$OUTDIR"; cd "$OUTDIR"
cat "$REFS" "$ISO" "$OG" > input.fasta

# 1. dereplicate full-length references
vsearch --derep_fulllength input.fasta --output deduped.fasta

# 2. align (adjustdirection also orients the outgroup)
mafft --auto --adjustdirection --thread 8 deduped.fasta > deduped_mafft.aln

# 3. ML tree (+ automatic rooting placeholder + SH-like support)
"$RAXML" -T 8 -m GTRGAMMA -s deduped_mafft.aln -n raxml -f d -p 12354 -\# "$NSEARCH"
"$RAXML" -T 8 -m GTRGAMMA -f I -t RAxML_bestTree.raxml -n root
"$RAXML" -T 8 -m GTRGAMMA -f J -p 12354 -t RAxML_rootedTree.root -n conf -s deduped_mafft.aln

# 4. add amplicon fragments to the reference alignment
mafft --auto --addfragments "$SHORT" --keeplength --thread 8 deduped_mafft.aln > addfragments.fasta

# 5. reference package + phylogenetic placement
taxit create -l 18S -P refpkg --aln-fasta deduped_mafft.aln \
  --tree-stats RAxML_info.raxml --stats-type RAxML --tree-file RAxML_fastTreeSH_Support.conf
pplacer -o placements.jplace -p --keep-at-most 20 -c refpkg addfragments.fasta

# 6. tables + trees for the figure scripts (unfiltered; figures filter in R)
guppy to_csv --point-mass --pp -o placements.csv placements.jplace
guppy tog -o tog.nwk placements.jplace
echo "DONE $OUTDIR  ->  tog.nwk, placements.csv"
