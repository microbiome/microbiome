data.dir <- system.file("extdata", package = "microbiome")

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  probedata <- as.matrix(tab)

  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  # Convert into phyloseq taxonomyTable format
  taxonomy <- tax_table(as.matrix(tab))     

require(affy)

# probesets
map <- unique(taxonomy[, c("species", "oligoID")]@.Data)
colnames(map) <- c("set", "probe")
rownames(map) <- map[, "probe"]

# Align the probe mappings with the data
# only use the probes that have mapping
coms <- intersect(rownames(probedata), rownames(map))
map <- map[coms,]
probedata <- probedata[coms,]

# Custom CDF
nsamples <- ncol(probedata)
nrow <- nrow(probedata)
ncol <- 1
sets <- unique(map[, "set"])
cdf <- lapply(sets, function(set) {
    tmp <- matrix(which(map[, "set"] == set), ncol = 1)
    colnames(tmp) <- "pm"
    return(tmp)
})
names(cdf) <- sets
cdfname <- list2env(cdf)

ab <- new("AffyBatch", exprs = probedata,
		 cdfName = "cdfname", 
            	 #phenoData = phenoData,
		 nrow = nrow,
		 ncol = ncol
            #annotation = cleancdfname(cdfname, addcdf = FALSE), 
            #protocolData = protocol,
	    #description = description, 
            #notes = notes
	    )

eset <- exprs(rma(ab, background = FALSE, normalize = FALSE))

# Not sure if this takes normalization out - check
eset <- exprs(rpa(ab, bg.method = "none", normalization.method = NULL))

## Check given species
#set <- sample(sets, 1)
## Summary
#x <- eset[set,]
## Probedata
#this.probes <- as.vector(map[map[, "set"] == set, "probe"])
#x2 <- probedata[this.probes,]
#cors <- cor(x, t(probedata), method = "spearman"); names(cors) <- rownames(probedata)
# Ranks wrt correct data
#sort(match(this.probes, rev(names(sort(cors)))))
# Ranks wrt random data
#sort(match(this.probes, rev(names(sample(cors)))))

# -----------------------------------------------


# Toydata example

nsamples <- 20
nrow <- 3
ncol <- 4

exprs <- matrix(runif(nrow * ncol * nsamples), nrow = nrow * ncol, ncol = nsamples)
colnames(exprs) <- paste("sample", 1:ncol(exprs), sep = "-")

# Custom CDF
maxIndx <- nrow * ncol
nsets <- 3
cdf <- lapply(seq(1, nsets), function(x) {
    tmp <- matrix(sample(maxIndx, 2), nrow = 2, ncol = 1)
    colnames(tmp) <- "pm"
    return(tmp)
})
names(cdf) <- paste("set", seq(1, nsets), sep = "")
cdf[["set1"]][, "pm"] <- c(2, 7)

v <- runif(nsamples)
exprs[2, ] <- v
exprs[7, ] <- v + runif(nsamples) * 0.01

cdfname <- list2env(cdf)

# OK !
ab <- new("AffyBatch", exprs = exprs,
		 cdfName = "cdfname", 
            	 #phenoData = phenoData,
		 nrow = nrow,
		 ncol = ncol
            #annotation = cleancdfname(cdfname, addcdf = FALSE), 
            #protocolData = protocol,
	    #description = description, 
            #notes = notes
	    )

indices2xy(30, nc = ncol)
xy2indices(x = 2, y = 3, nc = ncol)
exprs(ab)

# Summarize
# Jess !
exprs(rma(ab, background = FALSE, normalize = FALSE))

------------------------------------------

