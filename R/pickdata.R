pickdata <- function (x, otu.name) {

  taxa_are_rows <- NULL

  if (is.vector(x)) { 
    xxx <- x
  } else if (class(x) == "phyloseq") { 

    xx <- get_taxa(x)
    meta <- sample_data(x)
      
    if (!taxa_are_rows(x)) { xx <- t(xx)}

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