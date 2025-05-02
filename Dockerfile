# Use the Rocker image with R and Ubuntu Jammy
FROM rocker/r-ver:4.3.2 AS development

ARG STAFF=1000
ARG DEVELOPER=1000

# Image metadata
LABEL org.label-schema.license="MIT" \
      org.label-schema.vcs-url="https://gitlab.internal.sanger.ac.uk/sci/MAVEQC" \
      org.label-schema.vendor="MAVEQC-R Project" \
      maintainer="Cellular Informatics Scrum Team <cell-informatics@sanger.ac.uk>"

# Set working directory
WORKDIR /usr/src/MAVEQC-R

# Copy DESCRIPTION file with dependencies
COPY DESCRIPTION DESCRIPTION

# Create developer user, add it to staff group and provide permissions to working directory
# Install required dependencies
RUN groupadd -g ${STAFF} customstaff && \
    useradd -u ${DEVELOPER} -g customstaff -m -s /bin/bash developer \
    && apt-get update && apt-get install -y --no-install-recommends pandoc libssl-dev zlib1g-dev \
    && R -e "install.packages(c('BiocManager', 'devtools')); \
             devtools::install_deps('.'); \
             BiocManager::install(c('DESeq2', 'DEGreport', 'apeglm')); \
             devtools::install_version('MASS', version = '7.3-60', repos = 'https://cran.r-project.org'); \
             devtools::install_version('Matrix', version = '1.6-5', repos = 'https://cran.r-project.org'); \
             devtools::install_version('data.table', version = '1.14.10', repos = 'https://cran.r-project.org'); \
             install.packages('ggbeeswarm', repos = 'https://cran.r-project.org')" \
    && rm -rf /var/lib/apt/lists/*

# Initialise user
USER developer

# Default commandd
CMD ["/bin/bash"]
