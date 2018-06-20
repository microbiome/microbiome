library(microbiome)
#library(phyloseq)
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }

ressum <- read_hitchip("output", method = "sum")
resrpa <- read_hitchip("output", method = "rpa")
resfrpa <- read_hitchip("output", method = "frpa")

# Read precalculated
res.frpa <- read_profiling("output", method = "frpa")
res.rpa <- read_profiling("output", method = "rpa")
res.sum <- read_profiling("output", method = "sum")

# Check compatibility with old function
# frpa <- read.profiling(level = "L2", method = "frpa", data.dir = "output", log10 = TRUE, impute = TRUE) 

# Read and recalculate
res <- read_hitchip("output", method = "sum")

res1 <- res.sum$species
res2 <- res$abundance.table
res3 <- res.rpa$species
plot(log10(unlist(res1)), log10(unlist(res2)))
plot(log10(unlist(res1)), log10(unlist(res3)))


res1 <- res.rpa
res2 <- resrpa
plot(log10(unlist(res1$species)), log10(unlist(res2$abundance.table)))

