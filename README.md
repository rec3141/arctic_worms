# arctic_worms

Phylogenetic placement and biogeography of Arctic sea-ice flatworms (Acoela,
Rhabdocoela, Macrostomorpha) from 18S rRNA sequences. Cultured KBD isolates and
Sanger references are used to build rooted reference trees; V4 18S amplicon ASVs
are placed onto them, and isolate distributions are mapped across the Arctic.

## Repository layout

```
data/
  sequences/            reference, isolate, outgroup and amplicon 18S FASTAs
  trees/<group>/        placement outputs used by the figure scripts
                          tog.nwk        reference tree + placed ASVs (newick)
                          placements.csv per-placement stats (pendant length, LWR…)
  locations.tsv         sampling site GPS (one row per site)
  asv_site_summary.tsv  per-SITE mean relative abundance of each ASV *
  chimera_results.txt   uchime_ref output (chimera calls, Y/N)
scripts/
  run_pipeline.sh        derep → align → RAxML → pplacer placement (one group)
  prepare_site_summary.R aggregates raw DADA2 tables → data/asv_site_summary.tsv *
  plot_placement_tree.R  full rooted placement tree per group
  plot_kbd_panels.R      KBD overview (3 versions) + zoom panels per group
  plot_kbd_maps.R        per-isolate presence/absence + per-group pie & richness maps
figures/                 generated PDFs
```

\* Only **site-level** ASV summaries are distributed. The raw per-sample DADA2
tables are not included; `prepare_site_summary.R` documents how the summary was
produced from them.

## Reproduce

**Trees** (per group; needs `vsearch mafft raxmlHPC-PTHREADS-AVX taxit pplacer guppy gappa` on PATH):

```bash
scripts/run_pipeline.sh acoela_out \
  data/sequences/acoela_references.fst \
  data/sequences/acoela_outgroups.fasta \
  data/sequences/amplicons_platyhelminthes.fasta \
  data/sequences/acoela_isolates.fst 1000
# -> acoela_out/tog.nwk, placements.csv  (copy into data/trees/acoela/)
```
Rhabdocoela: swap in `rhabdocoela_*` + `catenulida_outgroups.fasta` + isolates.
Macrostomorpha: `macrostomorpha_references.fas` + `catenulida_outgroups.fasta`, no isolates.

Or rebuild **all three** groups at once and refresh `data/trees/`:

```bash
scripts/run_all.sh 1000        # nsearch; runs acoela, rhabdocoela, macrostomorpha
# in the container:
docker run --rm -v "$PWD":/work -w /work ghcr.io/rec3141/arctic_worms:latest scripts/run_all.sh 1000
```

**Figures** (R with `ape sf ggplot2 scatterpie rnaturalearth ggrepel`, plus
optional `showtext`+`sysfonts`; run from `scripts/`):

> `scripts/_fonts.R` renders tree-figure text as vector outlines when
> `showtext` and a sans font are available, so the PDFs display identically in
> every viewer (macOS Preview otherwise mangles cairo's subsetted fonts). It
> falls back to the device font if unavailable.


```bash
cd scripts
for g in acoela rhabdocoela macrostomorpha; do Rscript plot_placement_tree.R $g; done
for g in acoela rhabdocoela; do Rscript plot_kbd_panels.R $g; done
Rscript plot_kbd_maps.R
```

## Container / CI

The full tree-building toolchain (vsearch, mafft, RAxML, taxtastic, pplacer,
guppy, gappa) is packaged in a Docker image, built and smoke-tested by GitHub
Actions (`.github/workflows/docker.yml`) on every push and published to GHCR:

```bash
docker pull ghcr.io/rec3141/arctic_worms:latest

# reproduce a tree from the mounted repo (example: acoela, 1000 searches)
docker run --rm -v "$PWD":/work -w /work ghcr.io/rec3141/arctic_worms:latest \
  scripts/run_pipeline.sh acoela_out \
    data/sequences/acoela_references.fst \
    data/sequences/acoela_outgroups.fasta \
    data/sequences/amplicons_platyhelminthes.fasta \
    data/sequences/acoela_isolates.fst 1000
```

The image also carries the **R figure stack** (ape, sf, ggplot2, scatterpie,
rnaturalearth, ggrepel, showtext), so it reproduces trees *and* figures:

```bash
# full reproduction into the mounted repo
docker run --rm -v "$PWD":/work -w /work ghcr.io/rec3141/arctic_worms:latest bash -c \
  "scripts/run_all.sh 1000 && scripts/make_figures.sh"
```

Build it locally instead with `docker build -t arctic_worms .`. CI builds the
image and smoke-tests both the pipeline (macro, 2 searches → tree + placements)
and figure generation (R + sf/PROJ + showtext → PDFs).

## Methods notes

- **Outgroups:** Acoela rooted on Nemertodermatida (*Nemertoderma westbladi*
  AF327726, *Meara stichopi* AF119085); Rhabdocoela & Macrostomorpha on
  Catenulida (*Catenula lemnae* FJ196323, *Stenostomum leucops* FJ384832).
  The two per group form a paraphyletic grade, so trees are rooted on the
  branch between the outgroup and the ingroup.
- **ASV filtering (figures):** placements kept if best-placement
  `pendant_length < 0.2`; uchime-flagged chimeras removed. The strict gappa
  mass filter is intentionally *not* used, so confidently-placed but
  mass-split ASVs (near-identical reference clusters) are retained.
- **ASV ↔ isolate matching:** each ASV is assigned to its nearest KBD isolate
  if within 2% patristic distance (`ASV_MAX_DIV = 0.02`).
- ASV labels show the ASV number only (the `_Genus_species` in the FASTA
  headers is an automated DADA2/PR2 assignment, unreliable at species level).
