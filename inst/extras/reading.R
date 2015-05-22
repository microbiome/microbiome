library(microbiome)
#library(phyloseq)
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }

# Read precalculated
res.frpa <- read.profiling("output", method = "frpa")
res.rpa <- read.profiling("output", method = "rpa")
res.sum <- read.profiling("output", method = "sum")

# Read and recalculate
res <- read_hitchip("output", method = "sum")

res1 <- res$species
res2 <- res.sum$abundance.table
plot(log10(unlist(res1)), log10(unlist(res2)))

