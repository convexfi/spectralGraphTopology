FROM rocker/verse:latest
RUN Rscript -e "devtools::install_github('dppalomar/spectralGraphTopology')"
RUN Rscript -e "library(covr); codecov()"
