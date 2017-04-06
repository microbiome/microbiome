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
#' @param measures Default is ‘NULL’, meaning that all available
#'          alpha-diversity measures will be included. Alternatively, you
#'          can specify one or more measures as a character vector of
#'          measure names. Values include those supported in the
#'          phyloseq::estimate_richness function:
#'          "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher".
#'   In addition, the following measures are provided:
#'     "Evenness" (Pielou's index);
#'     "Dominance" (Number of species needed to cover 50% of the ecosystem);
#'     "Top" (Relative proportion of the most dominant species in [0,1]);
#'     "Rarity" (Relative proportion of the rare (non-core) species in [0,1]);
#'     "Core" (Relative proportion of the core species in [0,1]);
#'     "Gini" (Gini index).
#' @inheritParams core
#' @return A data.frame of samples x diversity indicators;
#'   except when split=FALSE, a vector of indices is returned.
#' @details This function returns the indices with the default choices for simplicity. See the individual functions for more options.
#' @examples 
#'   data(dietswap)
#'   d <- diversity(dietswap)
#' @export
#' @seealso rarity, core_abundance, dominance
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
diversity <- function(x, detection = 0, prevalence = 0, split = TRUE, measures = NULL) {

  res <- NULL
  measures1 <- setdiff(measures, c("Richness", "Evenness", "Dominance"))
  if (is.null(measures1) || length(measures1) > 0) {
    res <- estimate_richness(x, split = split, measures = measures1)
  } 

  if (("Evenness" %in% measures) || is.null(measures)) {
    # Shannon Diversity
    d <- estimate_richness(x, split = split, measures = "Shannon")$Shannon
    # normalize by Log richness to get Pielou's evenness
    r <- diversity(x, split = split, measures = "Observed")
    if (split) {r <- r$Observed} else {r <- unname(r)}
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

  if (("Dominance" %in% measures) || is.null(measures)) {

    do <- unname(dominance(x, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Dominance = do)
    } else {
      res$Dominance <- do
    }

  }


  if (("Gini" %in% measures) || is.null(measures)) {

    do <- unname(gini(x, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Gini = do)
    } else {
      res$Gini <- do
    }

  }

  # Calculate compositional abundance
  xcomp <- transform(x, "compositional")

  if (("Top" %in% measures) || is.null(measures)) {

    do <- unname(top_abundance(xcomp, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Top = do)
    } else {
      res$Top <- do
    }

  }

  if (("Rarity" %in% measures) || is.null(measures)) {

    do <- unname(rarity(xcomp, detection, prevalence, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Rarity = do)
    } else {
      res$Rarity <- do
    }

  }


  if (("Core" %in% measures) || is.null(measures)) {

    do <- unname(core_abundance(xcomp, detection, prevalence, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Core = do)
    } else {
      res$Core <- do
    }

  }

  # For total diversity, just return a vector of values
  # (as the table would have only 1 row; and giving a different
  #  output for this special case can potentially help avoid confusion)
  if (!split) {
    res <- unlist(res)
  }

  res

}



