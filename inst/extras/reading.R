library(microbiome)
#library(phyloseq)
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
res.frpa <- read.profiling("output", method = "frpa")
res.rpa <- read.profiling("output", method = "rpa")
res.sum <- read.profiling("output", method = "sum")

ressum <- read_hitchip("output", method = "sum")
resrpa <- read_hitchip("output", method = "rpa")
resfrpa <- read_hitchip("output", method = "frpa")

#plot(log10(unlist(res.rpa$species)), log10(unlist(res.frpa$species)))
#plot(log10(unlist(res.rpa$species)), log10(unlist(res.sum$species)))
plot(log10(unlist(res$species)), log10(unlist(res.sum)))

res1 <- res.rpa
res2 <- resrpa
plot(log10(unlist(res1$species)), log10(unlist(res2$abundance.table)))

