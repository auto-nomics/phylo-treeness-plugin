FROM docker.io/condaforge/miniforge3:latest@sha256:be01a3e63d4af20d8a9937d1e12fa299d1c572edf601b961180a59c78bd7250d

LABEL org.opencontainers.image.title="autonomics-phylo-treeness" \
      org.opencontainers.image.description="Pinned phylogenetics runtime: MAFFT alignment, ClipKIT trimming, IQ-TREE inference, PhyKIT treeness" \
      org.opencontainers.image.version="0.1.0" \
      org.opencontainers.image.source="https://github.com/auto-nomics/phylo-treeness-plugin" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# Exact bioconda pins; verify every executable resolves after the solve.
RUN mamba create -y -p /opt/phylo -c conda-forge -c bioconda \
        python=3.12 \
        mafft=7.526 \
        clipkit=2.14.0 \
        iqtree=3.1.4 \
        phykit=2.7.0 \
    && /opt/phylo/bin/mafft --version 2>&1 | grep -q "7.526" \
    && /opt/phylo/bin/clipkit --version | grep -q "2.14.0" \
    && /opt/phylo/bin/iqtree --version 2>&1 | grep -q "3.1.4" \
    && /opt/phylo/bin/phykit --help 2>&1 | grep -q "treeness" \
    && /opt/phylo/bin/python --version | grep -q "3.12"

ENV PATH=/opt/phylo/bin:${PATH}

RUN groupadd --gid 1001 autonomics \
    && useradd --uid 1001 --gid autonomics --no-create-home autonomics

ENV HOME=/tmp

USER 1001:1001

WORKDIR /work

ENTRYPOINT ["bash"]

CMD ["--version"]
