# Use the Bioconductor 3.18 image
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# 1. System dependencies for high-quality graphics and data handling
RUN apt-get update && apt-get install -y \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libpng-dev libjpeg-dev libtiff-dev libfreetype6-dev \
    libfribidi-dev libharfbuzz-dev libicu-dev \
    zlib1g-dev libfontconfig1-dev cmake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'clusterProfiler', \
    'enrichplot', \
    'GSVA', \
    'GSEABase', \
    'EnhancedVolcano', \
    'ComplexHeatmap', \
    'GOSemSim', \
    'treeio', \
    'ggtree' \
    ))"

# 3. Install CRAN packages (including dependencies for pairwise_termsim)
RUN R -e "install.packages(c( \
    'ggplot2', \
    'ggnewscale', \
    'circlize', \
    'stringr', \
    'data.table', \
    'magrittr', \
    'shadowtext', \
    'ggwordcloud' \
    ), repos='http://cran.rstudio.com/')"

WORKDIR /data
