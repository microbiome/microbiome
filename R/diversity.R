#' @title Estimate diversity
#' @description Diversity estimation. Augments the estimate_richness function of the phyloseq package
#' @param x \code{\link{phyloseq-class}} object 
#' @param split See help(phyloseq::estimate_richness)
#' @param measures See help(phyloseq::estimate_richness). In addition,
#'        the measure "Evenness" is provided (Pielou's index).
#' @param det.th detection threshold for observing taxa (absence / presence)
#' @return Vector containing relative proportions for each phylotype in 
#'         each sample 
#' @examples 
#'   #pseq <- download_microbiome("peerj32")$physeq
#'   #d <- estimate_diversity(pseq, det.th = 0)
#' @importFrom phyloseq estimate_richness
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
estimate_diversity <- function(x, split = TRUE, measures = NULL, det.th = 0) {

  res <- NULL
  measures1 <- setdiff(measures, c("Richness", "Evenness"))
  if (is.null(measures1) || length(measures1) > 0) {

    res <- estimate_richness(x, split = split, measures = measures1)

  } 

  if (("Evenness" %in% measures) || is.null(measures)) {
    # Shannon Diversity
    d <- estimate_richness(x, split = split, measures = "Shannon")$Shannon
    # normalize by Log richness to get Pielou's evenness
    r <- estimate_diversity(x, split = split, measures = "Observed")$Observed
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

    otu <- otu_table(x)@.Data

    # Calculate richness.
    # This simply indicates how many taxa are present in each sample
    # (exceed the detection threshold). This measure is sometimes used with
    # phylogenetic microarrays.
    ri <- colSums(otu > det.th)

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Observed = ri)
    } else {
      res$Observed <- ri
    }
  }

  res

}



