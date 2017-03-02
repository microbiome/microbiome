#' @title Estimate Diversity
#' @description Diversity estimation.
#'    Augments the estimate_richness function of the phyloseq package.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection detection for observing taxa (absence / presence).
#'   Used to determine the richness (Observed diversity) above this
#'   abundance threshold. Zero by default.
#' @param split (Optional). Logical. Should a separate set of richness
#'        estimates be performed for each sample? Or alternatively,
#'        pool all samples and estimate richness of the entire set.
#' @param measures (Optional). Default is ‘NULL’, meaning that all available
#'          alpha-diversity measures will be included. Alternatively, you
#'          can specify one or more measures as a character vector of
#'          measure names. Values must be among those supported in the
#'          phyloseq::estimate_richness function. These include
#'‘c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")’.
#' In addition, the measure "Evenness" is provided (Pielou's index).
#' @return A data.frame of samples x diversity indicators;
#'   except when split=FALSE, a vector of indices is returned.
#' @examples 
#'   data(dietswap)
#'   d <- diversity(dietswap)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
diversity <- function(x, detection = 0, split = TRUE, measures = NULL) {

  res <- NULL
  measures1 <- setdiff(measures, c("Richness", "Evenness"))
  if (is.null(measures1) || length(measures1) > 0) {
    res <- estimate_richness(x, split = split, measures = measures1)
  } 

  if (("Evenness" %in% measures) || is.null(measures)) {
    # Shannon Diversity
    d <- estimate_richness(x, split = split, measures = "Shannon")$Shannon
    # normalize by Log richness to get Pielou's evenness
    r <- diversity(x, split = split, measures = "Observed")$Observed
    e <- d/log(r)

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Evenness = e)
    } else {
      res$Evenness <- e
    }    
  }

  if (("Observed" %in% measures) || is.null(measures)) {

    # Pick the OTU data
    otu <- abundances(x)

    if (!split) {
      otu <- as.matrix(rowSums(otu), nrow = nrow(otu))
    }

    # Calculate richness.
    # This simply indicates how many taxa are present in each sample
    # (exceed the detection). This measure is sometimes used with
    # phylogenetic microarrays.
    ri <- colSums(otu > detection)

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Observed = ri)
    } else {
      res$Observed <- ri
    }

  }

  # For total diversity, just return a vector of values
  # (as the table would have only 1 row; and giving a different
  #  output for this special case can potentially help avoid confusion)
  if (!split) {
    unlist(res)
  }

  res

}



