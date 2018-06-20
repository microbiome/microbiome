library(microbiome)

# Old PET data, pre-calculated
d1 <- read.profiling(data.dir = "~/antagomir/pet14/data", level = "species", method = "frpa")

# New PET data, calculated on the fly
d2 <- read_hitchip(data.dir = "~/antagomir/pet14b/", method = "frpa")
d2 <- otu_table(d2$pseq)@.Data

cr <- intersect(rownames(d1), rownames(d2))
cc <- intersect(colnames(d1), colnames(d2))

plot(unlist(d1[cr, cc]), log10(unlist(d2[cr, cc])), pch = "."); abline(0,1)

cors <- c()
for (i in 1:nrow(frpa2)) {
  cors[[i]] <- cor(as.vector(log10(frpa2[cr, cc][i,])), as.vector(log10(frpa2b[cr, cc][i,])))
}
hist(cors)


