#fs <- list.files("~/HITChipDB/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }

dbuser <- "root"
dbpwd <- "fidipro"
dbname <- "phyloarray"
verbose = TRUE
host = '127.0.0.1'
port = 3307
summarization.methods = c("frpa", "sum")
which.projects = "AUTISM"

# Install dependencies
# git clone https://github.com/hadley/devtools.git

install.packages("sorvi")
install.packages("devtools")
source("http://www.bioconductor.org/biocLite.R")
biocLite("WGCNA")

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



