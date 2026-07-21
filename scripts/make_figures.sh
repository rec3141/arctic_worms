#!/bin/bash
# Regenerate every figure from data/trees/ + data/ into figures/.
# Requires R with: ape sf ggplot2 scatterpie rnaturalearth ggrepel showtext.
set -euo pipefail
cd "$(dirname "$0")"                 # scripts/
mkdir -p ../figures

for g in acoela rhabdocoela macrostomorpha; do Rscript plot_placement_tree.R "$g"; done
for g in acoela rhabdocoela; do Rscript plot_kbd_panels.R "$g"; done
Rscript plot_kbd_maps.R

echo "DONE — figures written to figures/"
