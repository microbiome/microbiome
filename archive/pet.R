library(microbiome)

# New pipeline
frpa1 <- read.profiling(data.dir = "PET", level = "species", method = "frpa")
frpa2 <- read.profiling(data.dir = "PET", level = "species", method = "sum")

# New pipeline 2
sum.new2 <- read_hitchip("PET", method = "sum")
frpa.new2 <- read_hitchip("PET", method = "frpa") 
rpa.new2 <- read_hitchip("PET", method = "rpa") 

# Old pipeline
frpa.old <- read.profiling("~/tmp/pet14/data", method = "frpa")
sum.old <- read.profiling("~/tmp/pet14/data", method = "sum")

# Old and New match perfect at species level
level <- "species"
s <- colnames(sum.old[[level]])
taxa <- intersect(rownames(sum.new[[level]]), rownames(sum.old[[level]]))
plot(log10(unlist(sum.new[[level]][taxa,s])), log10(unlist(sum.old[[level]][taxa,s])))
plot(log10(unlist(frpa.new[[level]][taxa,s])), log10(unlist(frpa.old[[level]][taxa,s])))

level <- "species"
par(mfrow = c(2,2))
plot(log10(unlist(sum.new[[level]])), log10(unlist(sum.new2$abundance.table)))
plot(log10(unlist(frpa.new[[level]])), log10(unlist(frpa.new2$abundance.table)))
plot(log10(unlist(frpa.new[[level]])), log10(unlist(sum.new2$abundance.table)))

# FIXME: make this work also with old output
# Read and recalculate
res <- read_hitchip("PET", method = "sum")
pseq <- res$pseq
otu <- res$abundance.table

#  pseq@otu_table
dim(tax_glom(pseq, "L2")@otu_table)
dim(tax_glom(pseq, "L1")@otu_table)

plot(unlist(log10(tax_glom(pseq, "species")@otu_table)), unlist(log10(pseq@otu_table)))
sum.new <- read.profiling("PET", method = "sum")


#---------------------------------------------------------------

sum.old <- read.profiling("~/tmp/pet14/data", method = "sum")
sum.new <- read.profiling("PET", method = "sum")
pseq0 <- read_hitchip("PET", method = "sum", detection.threshold = 0)$pseq
pseq <- read_hitchip("PET", method = "sum", detection.threshold = 10^1.8)$pseq
pseq.L20 <- aggregate_taxa(pseq0, level = "L2")
pseq.L2 <- aggregate_taxa(pseq, level = "L2")

taxa <- intersect(rownames(L2), rownames(sum.new$L2))

coms <- colnames(sum.old$L2[taxa,])

old <- log10(sum.old$L2[taxa,coms])
new.pre <- log10(sum.new$L2[taxa,coms])
new0 <- log10(otu_table(pseq.L20)@.Data[taxa,coms])
new.th <- log10(otu_table(pseq.L2)@.Data[taxa,coms])

par(mfrow = c(2,2))
plot(unlist(old[, coms]), unlist(new0[, coms]))
plot(unlist(old[, coms]), unlist(new.th[, coms]))
plot(unlist(new.th[, coms]), unlist(new.pre[, coms]))







