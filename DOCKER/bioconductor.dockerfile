# Use the Bioconductor 3.18 image
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# ------------------------------------------------------------------------------
# LAYER 1: System Dependencies (apt-get)
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libpng-dev libjpeg-dev libtiff-dev libfreetype6-dev \
    libfribidi-dev libharfbuzz-dev libicu-dev \
    zlib1g-dev libfontconfig1-dev cmake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# LAYER 2: Human Databases
# ------------------------------------------------------------------------------
RUN R -e "BiocManager::install('org.Hs.eg.db')"

# ------------------------------------------------------------------------------
# LAYER 3: Rat Databases
# ------------------------------------------------------------------------------
RUN R -e "BiocManager::install('org.Rn.eg.db')"

# ------------------------------------------------------------------------------
# LAYER 4: Core Computational Engines & Orthology
# ------------------------------------------------------------------------------
RUN R -e "BiocManager::install(c( \
    'DESeq2', \
    'clusterProfiler', \
    'GSVA', \
    'GSEABase', \
    'EnhancedVolcano', \
    'biomaRt' \
    ))"

# ------------------------------------------------------------------------------
# LAYER 5: Foundation CRAN Packages
# ------------------------------------------------------------------------------
RUN R -e "install.packages(c( \
    'data.table', \
    'dplyr', \
    'stringr', \
    'magrittr', \
    'ggplot2', \
    'ggnewscale', \
    'msigdbr' \
    ), repos='http://cran.rstudio.com/')"

# ------------------------------------------------------------------------------
# LAYER 6: Complex Visualization Libraries
# ------------------------------------------------------------------------------
RUN R -e "BiocManager::install(c( \
    'ComplexHeatmap', \
    'enrichplot', \
    'GOSemSim', \
    'treeio', \
    'ggtree' \
    ))"

RUN R -e "install.packages(c( \
    'circlize', \
    'shadowtext', \
    'ggwordcloud', \
    'ggupset', \
    'pheatmap', \
    'RColorBrewer' \
    ), repos='http://cran.rstudio.com/')"

# ------------------------------------------------------------------------------
# LAYER 7: Volatile Layer
# ------------------------------------------------------------------------------
RUN R -e "install.packages(c( \
    'ape', \
    'ggrepel' \
    ), repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install('EnsDb.Hsapiens.v86')"
RUN R -e "BiocManager::install('EnsDb.Rnorvegicus.v79')"

# ------------------------------------------------------------------------------
# LAYER 8: Final Additions
# ------------------------------------------------------------------------------
RUN R -e "install.packages(c( \
    'scatterplot3d', \
    'tidyr' \
    ), repos='http://cran.rstudio.com/')"

# ------------------------------------------------------------------------------
# LAYER 9: PPI Network
# ------------------------------------------------------------------------------
RUN R -e "install.packages(c( \
    'R.utils', \
    'ggraph', \
    'igraph' \
    ), repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install('DOSE')"

# ------------------------------------------------------------------------------
# LAYER 10: Excel Generation (NEW)
# ------------------------------------------------------------------------------
RUN R -e "install.packages('openxlsx', repos='http://cran.rstudio.com/')"

# ------------------------------------------------------------------------------
# LAYER 11: Statistical Reporting (Global Subtypes)
# ------------------------------------------------------------------------------
RUN R -e "install.packages(c( \
    'vegan', \
    'knitr' \
    ), repos='http://cran.rstudio.com/')"

# ------------------------------------------------------------------------------
# LAYER 12: High-Sensitivity Stats (Car + Limma)
# ------------------------------------------------------------------------------
RUN R -e "install.packages('car', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('limma')"

WORKDIR /data
