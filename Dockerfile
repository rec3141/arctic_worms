# Reproducible environment for the arctic_worms RAxML / pplacer pipeline.
# All tools come from conda-forge / bioconda (+ taxtastic from pip).
# Built for linux/amd64 (bioconda pplacer is linux-64 only).
FROM --platform=linux/amd64 condaforge/miniforge3:latest

# --- pipeline tools -------------------------------------------------------
RUN mamba install -y -n base -c conda-forge -c bioconda \
        "python=3.10" vsearch mafft raxml pplacer gappa \
 && pip install --no-cache-dir taxtastic \
 && mamba clean -a -y

WORKDIR /arctic_worms
COPY . /arctic_worms
RUN chmod +x scripts/*.sh scripts/*.R

# --- sanity check: every tool resolves ------------------------------------
RUN vsearch --version 2>&1 | head -1 \
 && mafft --version 2>&1 | head -1 \
 && (ls /opt/conda/bin | grep -i 'raxmlHPC' | head) \
 && pplacer --version \
 && guppy --version \
 && gappa --version \
 && taxit --help >/dev/null \
 && echo "== all pipeline tools present =="

CMD ["bash"]
