FROM bioconductor/bioconductor_docker:RELEASE_3_19

MAINTAINER pl.germain@gmail.com

WORKDIR /home/build/package

COPY . /home/build/package 

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN pip3 install snakemake

RUN Rscript -e 'BiocManager::install(c("BiocParallel", "ComplexHeatmap", "cowplot", "data.table", "edgeR", "ggplot2", "ggrastr", "IHW", "muscat", "patchwork", "RColorBrewer", "SingleCellExperiment", "SummarizedExperiment", "RColorBrewer", "IHW", "scuttle", "SEtools", "circlize", "viridis", "dplyr", "reshape2", "ggrastr", "matrixStats", "scater", "ComplexHeatmap", "scales"));'
RUN Rscript -e 'BiocManager::install("PRROC")'
RUN Rscript -e 'BiocManager::install("HelenaLC/muscat", ref="ff6d5a020c2171deaae67c9686a7974e296ca78d", build_vignettes=FALSE)'

