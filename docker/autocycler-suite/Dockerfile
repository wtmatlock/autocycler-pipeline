FROM --platform=linux/amd64 mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

USER root 

RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    python==3.8.20 \
    autocycler==0.4.0 \
    canu==2.3 \
    flye==2.9.5 \
    raven-assembler==1.8.3 \
    miniasm==0.3 \
    minipolish==0.2.0 \
    minimap2==2.29 \
    racon==1.5.0 \
    seqtk==1.5 \
    any2fasta==0.4.2 \
    plassembler==1.7.0 \
    unicycler==0.5.0 \
    file==5.46 \
    && micromamba clean --all --yes
 
RUN apt-get update && apt-get install -y --no-install-recommends \
    gawk \
    grep \
    procps \
    sed \
    coreutils \
    util-linux \
    bsdmainutils \
    && rm -rf /var/lib/apt/lists/*

# Make sure every layer runs with the activated base env
SHELL ["micromamba", "run", "-n", "base", "/bin/bash", "-lc"]

# Default command to autocycler help
CMD ["micromamba", "run", "-n", "base", "autocycler", "--help"]