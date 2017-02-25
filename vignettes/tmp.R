# Get example data
library(microbiome)
data("peerj32")
pseq <- transform_phyloseq(peerj32$phyloseq, "hellinger")
otu <- taxa_abundances(transform_phyloseq(pseq, "log10"))
meta <- pseq_metadata(pseq)
meta$group <- meta[["gender"]]


### PERMANOVA

# Use relative abundances for simpler visualizations
# PERMANOVA: samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ group, data=meta, permutations=99, method = "euclidean")
pv.permanova <- as.data.frame(permanova$aov.tab)["group", "Pr(>F)"]
efs.permanova <- permanova$coefficients["group1",] # Effect sizes

# ----------------------------

# Compare the two groups with limma
library(limma)

# Prepare the design matrix which states the groups for each sample
# in the otu
design <- cbind(intercept = 1, Grp2vs1 = meta$group)
rownames(design) <- rownames(meta)
design <- design[colnames(otu), ]

# NOTE: results and p-values are given for all groupings in the design matrix
# Now focus on the second grouping ie. pairwise comparison
coef.index <- 2
     
# Fit the limma model
fit <- lmFit(otu, design)
fit <- eBayes(fit)
pv.limma <- fit$p.value[, 2]
efs.limma <-  fit$coefficients[, "Grp2vs1"]

# -------------------------------------------------

# Compare the two groups with t-test
library(dplyr)
male.samples <- suppressWarnings(dplyr::filter(meta, gender == "male")$sample)
female.samples <- suppressWarnings(dplyr::filter(meta, gender == "female")$sample)

# Compare p-values between limma and permanova
taxa <- rownames(otu)
#plot(efs.limma[taxa], efs.permanova[taxa])
plot(log10(1 + pv.limma[taxa]), efs.permanova[taxa])
#abline(0,1,lty = 2, xlim = c(0,1), ylim = c(0,1))


