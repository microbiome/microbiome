# Use this script to generate the example data sets
# for the package

source("download.R")

# Download the required R packages and then the HITChip Atlas data set
library("rdryad")
library("microbiome")
atlas1006 <- download_microbiome("atlas1006", "phyloseq")
save(atlas1006, file = "../../data/atlas1006.rda")

# -------------------------------------------------------

library(microbiome)
dietswap <- download_microbiome("dietswap")
save(dietswap, file = "../../data/dietswap.rda")

# ----------------------------------------------

library(microbiome)
peerj32 <- download_microbiome("peerj32")
save(peerj32, file = "../../data/peerj32.rda")
