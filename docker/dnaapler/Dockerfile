FROM --platform=linux/amd64 mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

# Install dnaapler and dependencies into the base environment
RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    python=3.10 \
    dnaapler \
    mmseqs2==13.45111 \
    just \
    poetry \
    ripgrep \
    pyrodigal \
    && micromamba clean --all --yes

# Make sure every layer runs with the activated base env
SHELL ["micromamba", "run", "-n", "base", "/bin/bash", "-lc"]

# Default command to dnaapler help
CMD ["micromamba", "run", "-n", "base", "dnaapler", "--help"]