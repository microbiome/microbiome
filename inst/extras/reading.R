library(microbiome)
#library(phyloseq)
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
res.frpa <- read.profiling("temp2", method = "frpa")
res.rpa <- read.profiling("temp2", method = "rpa")
res.sum <- read.profiling("temp2", method = "sum")
#res.sum <- read_hitchip("temp2", method = "sum")
#plot(log10(unlist(res.rpa$species)), log10(unlist(res.frpa$species)))
#plot(log10(unlist(res.rpa$species)), log10(unlist(res.sum$species)))
plot(log10(unlist(res$species)), log10(unlist(abu)))

