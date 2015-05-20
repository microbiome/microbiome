data.dir <- system.file("extdata", package = "microbiome")

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  probedata <- tab

  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  # Convert into phyloseq taxonomyTable format
  taxonomy <- tax_table(as.matrix(tab))     

# ------------------------------------------

# probedata for a single sample
p <- probedata[,1]
names(p) <- rownames(probedata)

# probesets
map <- unique(taxonomy[, c("species", "oligoID")]@.Data)
rownames(map) <- NULL

# create probesets
sets <- split(map[, "oligoID"], map[, "species"])
allSetDat <- lapply(sets, function (probes) {m <- as.matrix(p[probes]); colnames(m) <- "pm"; m})

require(affy)
hgu2 <- list2env(allSetDat)
celDat@cdfName <- "hgu2"
rma4 <- exprs(rma(celDat))

# -----------------------------------------------

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

