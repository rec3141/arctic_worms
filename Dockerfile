# Reproducible environment for the arctic_worms pipeline AND figures.
# Pipeline tools (vsearch/mafft/RAxML/pplacer/gappa/taxtastic) + R figure stack
# (ape/sf/ggplot2/scatterpie/rnaturalearth/showtext) from conda-forge/bioconda.
# Built for linux/amd64 (bioconda pplacer is linux-64 only).
FROM --platform=linux/amd64 condaforge/miniforge3:latest

# --- pipeline tools + R figure stack (one env) ----------------------------
RUN mamba install -y -n base -c conda-forge -c bioconda \
        "python=3.10" vsearch mafft raxml pplacer gappa \
        r-base r-ape r-ggplot2 r-sf r-scatterpie r-ggrepel \
        r-showtext r-sysfonts r-rnaturalearth r-rnaturalearthdata \
 && pip install --no-cache-dir taxtastic \
 && mamba clean -a -y

# --- fonts so showtext can outline text -----------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        fonts-liberation fonts-dejavu-core \
 && rm -rf /var/lib/apt/lists/*

# GDAL/PROJ data paths (conda activation scripts aren't sourced under Rscript)
ENV PROJ_DATA=/opt/conda/share/proj \
    PROJ_LIB=/opt/conda/share/proj \
    GDAL_DATA=/opt/conda/share/gdal

WORKDIR /arctic_worms
COPY . /arctic_worms
RUN chmod +x scripts/*.sh scripts/*.R && mkdir -p figures   # figures/ is .dockerignore'd

# --- sanity: pipeline tools + R figure stack all resolve ------------------
RUN vsearch --version 2>&1 | head -1 \
 && (ls /opt/conda/bin | grep -i 'raxmlHPC' | head -1) \
 && pplacer --version && guppy --version && gappa --version && taxit --help >/dev/null \
 && Rscript -e 'suppressMessages({library(ape);library(sf);library(ggplot2);library(scatterpie);library(rnaturalearth);library(ggrepel);library(showtext)}); cat("R figure stack OK\n")'

CMD ["bash"]
