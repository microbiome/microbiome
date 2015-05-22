library(microbiome)
#library(phyloseq)
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
res.frpa <- read.profiling("output", method = "frpa")
res.rpa <- read.profiling("output", method = "rpa")
res.sum <- read.profiling("output", method = "sum")

res <- read_hitchip("output", method = "sum")

#plot(log10(unlist(res.rpa$species)), log10(unlist(res.frpa$species)))
#plot(log10(unlist(res.rpa$species)), log10(unlist(res.sum$species)))
plot(log10(unlist(res$species)), log10(unlist(res.sum)))

