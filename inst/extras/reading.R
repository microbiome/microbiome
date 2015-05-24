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

res1 <- res.sum$species
res2 <- res$abundance.table
res3 <- res.rpa$species
plot(log10(unlist(res1)), log10(unlist(res2)))
plot(log10(unlist(res1)), log10(unlist(res3)))

