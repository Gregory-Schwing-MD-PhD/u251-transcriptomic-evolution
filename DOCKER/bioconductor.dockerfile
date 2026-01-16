# Use the Bioconductor 3.18 image
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# 1. System dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libpng-dev libjpeg-dev libtiff-dev libfreetype6-dev \
    libfribidi-dev libharfbuzz-dev libicu-dev \
    zlib1g-dev libfontconfig1-dev cmake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Install Bioconductor packages
# Added DESeq2 explicitly just in case
RUN R -e "BiocManager::install(c( \
    'DESeq2', \
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

# 3. Install CRAN packages
# ADDED: ape, ggrepel, pheatmap, RColorBrewer, dplyr
RUN R -e "install.packages(c( \
    'ggplot2', \
    'ggnewscale', \
    'circlize', \
    'stringr', \
    'data.table', \
    'magrittr', \
    'shadowtext', \
    'ggwordcloud', \
    'ggupset', \
    'ape', \
    'ggrepel', \
    'pheatmap', \
    'RColorBrewer', \
    'dplyr' \
    ), repos='http://cran.rstudio.com/')"

WORKDIR /data
