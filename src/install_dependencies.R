## R package dependencies -------------------------------------------------------
install.packages(c("matrixStats",
                   "ggplot2",
                   "MASS",
                   "scales",
                   "devtools",
                   "VennDiagram",
                   "gridExtra")
                 repos = 'https://watson.nci.nih.gov/cran_mirror/')

## Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("GeneOverlap")

## SeqRetriever installation
devtools::install_github("hilldr/SeqRetriever/SeqRetriever")
