# medaka from bioconda with pyabpoa pinned to the version medaka checks for at run time. The bioconda recipes of
# medaka 2.1.1 do not pin pyabpoa, so every published 2.1.1 image carries a newer one and `medaka tandem` refuses
# to start ("pyabpoa == 1.5.1 required"). The resolved environment is exported next to the tool.
# Debian 12 (glibc 2.36): newer glibc refuses the PyTorch build medaka 2.1.1 resolves to ("cannot enable executable stack")
FROM mambaorg/micromamba:2.3.3-debian12-slim
ARG MEDAKA_VERSION
ARG PYABPOA_VERSION
ARG PYTHON_VERSION
RUN micromamba install -y -n base -c conda-forge -c bioconda --channel-priority strict \
        "medaka=${MEDAKA_VERSION}" "pyabpoa=${PYABPOA_VERSION}" "python=${PYTHON_VERSION}" procps-ng \
    && micromamba env export -n base --explicit > /opt/conda/medaka-${MEDAKA_VERSION}.lock.txt \
    && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN medaka --version && python -c "import pyabpoa, importlib.metadata as m; print('pyabpoa', m.version('pyabpoa'))"
ENV PATH=/opt/conda/bin:$PATH
