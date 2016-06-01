pickdata <- function (x, otu.name) {

  if (!class("a") == "character") {
    stop("Provide proper variable name for pickdata function.")
  }

  if (is.vector(x)) {
  
    xxx <- x
    
  } else if (class(x) == "phyloseq") { 

    xx <- taxa_abundances(x)
    meta <- sample_data(x)

    # If OTU name not in otu data then try metadata
    if (otu.name %in% rownames(xx)) {
      xxx <- as.vector(xx[otu.name,])
    } else if (otu.name %in% colnames(meta)) {
      xxx <- unlist(meta[, otu.name])
    } else {   
      stop(paste("The provided variable name", otu.name, "cannot be found from OTUs or metadata"))
    }
  }
  xxx
}