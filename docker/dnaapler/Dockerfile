FROM --platform=linux/amd64 mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

USER root 

RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    python==3.10.14 \
    dnaapler==1.2.0 \
    mmseqs2==13.45111 \
    just==1.40.0 \
    poetry==2.1.3 \
    ripgrep==14.1.1 \
    pyrodigal==3.6.3 \
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

# Default command to dnaapler help
CMD ["micromamba", "run", "-n", "base", "dnaapler", "--help"]