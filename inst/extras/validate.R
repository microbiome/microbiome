library(microbiome)
library(HITChipDB)

# Read probe level data
data.dir <- system.file("extdata", package = "microbiome")
res <- read_hitchip(data.dir, method = "sum", detection.threshold = 0)
probedata <- res$probedata
taxonomy <- res$taxonomy

# Preprocess with microbiome
level <- "species"
method <- "sum"
downweight.ambiguous.probes <- TRUE
verbose <- TRUE
data.dir <- NULL
probe.parameters = NULL

dm <- summarize_probedata(probedata = probedata,
      	 	             taxonomy = taxonomy,
      	 		     level = level,
			     method = method)

# ---------------------------------------------------------

# Preprocess with HITChipDB
dh <- summarize.sum(taxonomy, level, probedata, downweight.ambiguous.probes = TRUE)$abundance.table

# ----------------------------------------------------

# Compare
comsr <- intersect(rownames(dh), rownames(dm))
comsc <- intersect(colnames(dh), colnames(dm))
dhs <- dh[comsr, comsc]
dms <- dm[comsr, comsc]
plot(log10(as.vector(dhs)), log10(as.vector(dms))); abline(0,1)

# Ambiguous species giving a different outcome
amb <- names(which(sort(rowSums(abs(dhs - dms) > 0)) > 0))
set <- amb[[1]]