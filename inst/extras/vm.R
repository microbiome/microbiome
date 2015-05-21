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
biocLite("WGCNA")
biocLite("BiocStyle")
biocLite("knitr")
biocLite("rmarkdown")
biocLite("phyloseq")
biocLite("RPA")

library(devtools)
install_github("microbiome/microbiome")
install_github("microbiome/HITChipDB")

# -------------------------------------------

library(HITChipDB)
dbuser <- "root"
dbpwd <- "fidipro"
dbname <- "phyloarray"
verbose = TRUE
host = '127.0.0.1'
port = 3307
summarization.methods = c("frpa", "sum", "rpa")
which.projects = "AUTISM"

#fs <- list.files("~/HITChipDB/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
res <- run.profiling.script(dbuser, dbpwd, dbname, verbose = FALSE, host = host, port = port, summarization.methods = summarization.methods, which.projects = which.projects)

# -------------------------------

library(microbiome)
#library(phyloseq)
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
res <- read_hitchip("~/R/temp2", method = "frpa")




