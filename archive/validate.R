#' @title Validate Phyloseq
#' @description Validate phyloseq object.
#' @param x phyloseq object
#' @details Checks that the abundances and sample_data have exactly same samples.
#' @return A validated and polished phyloseq object
#' @export
#' @examples
#'   data(dietswap)
#'   validate(dietswap)
validate <- function (x) {

  validated <- TRUE		  

  dat <- t(abundances(x))
  meta <- sample_data(x)		         
  coms <- intersect(rownames(dat), rownames(meta))

  if (length(coms) < 2) {
    validated <- FALSE
    warning("Check that the abundances and sample_data have more than 1 samples 
             in common")
    return(validated)
  } else {
    # Include only the common samples    
    x@sam_data <- meta[coms,]
    otu_table(x) <- otu_table(t(dat[coms,]), taxa_are_rows = TRUE)
  }

  # Return validated and polished phyloseq object
  if (validated) {
    x
  }
}
