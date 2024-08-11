FROM r-base:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
libcurl4-openssl-dev \
libxml2-dev \
libssl-dev \
libhdf5-dev \
libharfbuzz-dev \
libfribidi-dev \
libudunits2-dev \
libgdal-dev \
libproj-dev \
libgeos-dev \
libfontconfig1-dev \
&& rm -rf /var/lib/apt/lists/*

#RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('rjson')"
RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(version = "3.18")'
# Installing Bioconductor core packages:
#RUN Rscript -e 'BiocManager::install(c("GenomicFeatures", "AnnotationDbi"), dependencies = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install(c("GenomicFeatures"), dependencies = TRUE, ask = FALSE)'
RUN Rscript -e 'install.packages("argparse")'