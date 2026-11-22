# medaka from bioconda with a chosen pyabpoa. The bioconda recipes of medaka 2.1.1 do not pin pyabpoa, so every
# published 2.1.1 image carries a newer one and `medaka tandem` refuses to start ("pyabpoa == 1.5.1 required").
# - PYABPOA_VERSION=1.5.1: the version medaka 2.1.1 checks for. On HG002 this combination fails silently in the
#   consensus stage (nanoporetech/medaka#583: an abpoa problem, fixed with abpoa 1.5.5).
# - PYABPOA_VERSION=1.5.5 (the working build): the one change to medaka's code is the version check in
#   medaka/tandem/tandem.py (`abpoa_required_version = '1.5.1'` -> the chosen version), the same change the reporters of
#   #583 made. Nothing else is touched. Listed as A-medaka-1 in docs/ADAPTATIONS.md.
# The resolved environment is exported next to the tool.
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
RUN if [ "${PYABPOA_VERSION}" != "1.5.1" ]; then \
        f=$(python -c "import medaka.tandem.tandem as t; print(t.__file__)") \
        && test "$(grep -c "abpoa_required_version = '1.5.1'" "$f")" = 1 \
        && sed -i "s/abpoa_required_version = '1.5.1'/abpoa_required_version = '${PYABPOA_VERSION}'/" "$f" \
        && grep -n "abpoa_required_version = " "$f"; \
    fi
RUN medaka --version && python -c "import pyabpoa, importlib.metadata as m; print('pyabpoa', m.version('pyabpoa'))"
ENV PATH=/opt/conda/bin:$PATH
