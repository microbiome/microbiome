#fs <- list.files("~/HITChipDB/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }

# Mahdollisesti riittaa:
# sudo apt-get update
# sudo apt-get install libxml2-dev

# sudo apt-get upgrade r-base
# wget http://cran-mirror.cs.uu.nl/src/base/R-3/R-3.2.0.tar.gz
# tar -zxvf R-3.2.0.tar.gz
# ./configure
# make
# sudo make install
# install.packages("devtools")

source("http://www.bioconductor.org/biocLite.R")
biocLite("RPA")
biocLite("phyloseq")
biocLite("WGCNA")
biocLite("BiocStyle")
biocLite("knitr")
biocLite("rmarkdown")

library(devtools)
install_github("microbiome/microbiome")
install_github("microbiome/HITChipDB")

# -------------------------------------------

dbuser <- "root"
dbpwd <- "fidipro"
dbname <- "phyloarray"
verbose = TRUE
host = '127.0.0.1'
port = 3307
summarization.methods = c("frpa", "sum")
which.projects = "AUTISM"

install.packages("sorvi")
library(HITChipDB)
library(RMySQL)
library(DBI)
library(tcltk)
library(phyloseq)
library(RPA)
fs <- list.files("~/HITChipDB/R/", full.names = TRUE)
for (f in fs)  { source(f) }
res <- run.profiling.script(dbuser, dbpwd, dbname, verbose = TRUE, host = host, port = port, summarization.methods = c("frpa", "sum"), which.projects = "AUTISM")

# -------------------------------

library(phyloseq)
fs <- list.files("~/microbiome/R/", full.names = TRUE)
for (f in fs)  { source(f) }
res <- read.profiling("~/R/temp2")



