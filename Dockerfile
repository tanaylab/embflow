FROM bioconductor/bioconductor_docker:RELEASE_3_12

# Install rpm dependencies
RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev libgtk2.0-dev libcairo2-dev libxt-dev xvfb xauth xfonts-base vim && rm -rf /var/lib/apt/lists/*


RUN R -e 'remotes::install_github("tanaylab/tgconfig")'
RUN R -e 'remotes::install_github("tanaylab/tglkmeans")'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("pheatmap")'
RUN R -e 'install.packages("gridExtra")'
RUN R -e 'install.packages("umap")'
RUN R -e 'install.packages("Matrix")'
RUN R -e 'install.packages("shape")'
RUN R -e 'install.packages("qlcMatrix")'
RUN R -e 'install.packages("ggrepel")'
RUN R -e 'BiocManager::install("tanaylab/metacell")'
RUN R -e 'BiocManager::install("lpsymphony")'

RUN git clone https://github.com/tanaylab/embflow.git

WORKDIR /embflow

RUN R -e 'source("scripts/download_data.r")'

# Run R
CMD ["R"]

