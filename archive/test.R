

library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")

method <- "sum"

so <- read.csv(paste("~/Rpackages/microbiome/microbiome/inst/extdata/species-", method, ".tab", sep = ""), sep = "\t", row.names = 1)

pseq <- read_hitchip(data.directory, method = method)$pseq
sn <- otu_table(pseq)@.Data

cro <- intersect(rownames(sn), rownames(so))
cco <- intersect(colnames(sn), colnames(so))
plot(log10(unlist(so[cro, cco])), log10(unlist(sn[cro, cco])))

cors <- c()
for (i in 1:length(cro)) {
  tax <- cro[[i]]
  cors[[tax]] <- cor(log10(unlist(so[tax, ])), log10(unlist(sn[tax, ])))
}
hist(cors)

# --------------------------------

# The poorest correlate
tax <- names(sort(cors))[[1]]


