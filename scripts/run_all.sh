#!/bin/bash
# ---------------------------------------------------------------------------
# Re-run the full RAxML/pplacer placement pipeline for all three taxon groups
# and refresh the trees used by the figure scripts (data/trees/<group>/).
#
# In the container (mount the repo so outputs persist to the host):
#   docker run --rm -v "$PWD":/work -w /work ghcr.io/rec3141/arctic_worms:latest \
#     scripts/run_all.sh [nsearch]        # nsearch default 1000
# ---------------------------------------------------------------------------
set -euo pipefail
NSEARCH=${1:-1000}
HERE=$(cd "$(dirname "$0")/.." && pwd); cd "$HERE"      # repo root
SEQ="$HERE/data/sequences"; AMP="$SEQ/amplicons_platyhelminthes.fasta"

run() {  # <group> <references> <outgroups> <isolates|/dev/null>
  local g=$1
  echo "======== $g  (nsearch=$NSEARCH) ========"
  rm -rf "out_$g"
  scripts/run_pipeline.sh "out_$g" "$SEQ/$2" "$SEQ/$3" "$AMP" "$4" "$NSEARCH"
  mkdir -p "data/trees/$g"
  cp "out_$g/tog.nwk" "out_$g/placements.csv" "data/trees/$g/"
  echo ">>> refreshed data/trees/$g/{tog.nwk,placements.csv}"
}

run acoela         acoela_references.fst         acoela_outgroups.fasta     "$SEQ/acoela_isolates.fst"
run rhabdocoela    rhabdocoela_references.fst    catenulida_outgroups.fasta "$SEQ/rhabdocoela_isolates.fst"
run macrostomorpha macrostomorpha_references.fas catenulida_outgroups.fasta /dev/null

echo "ALL DONE — regenerate figures next:  cd scripts && Rscript plot_*.R"
