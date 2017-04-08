#' @title Global Ecosystem State Variables 
#' @description Global indicators of the ecoystem state, including richness, evenness, diversity, and other indicators
#' @param x \code{\link{phyloseq-class}} object
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
#'     "Richness" (number of unique taxa that give non-zero signal); 
#'     "Evenness" (Pielou's index);
#'     "Dominance" (Number of species needed to cover 50% of the ecosystem);
#'     "Top_Abundance" (Relative proportion of the most dominant species in [0,1]);
#'     "Rarity" (Relative proportion of the rare (non-core) species in [0,1]) - this complement (1-x) of the Core_Abundance
#'     "Low_Abundance" (Relative proportion of the least abundant species, below the detection level of 0.2%); 
#'     "Core_Abundance" (Relative proportion of the core species that exceed detection level 0.2% in over 50% of the samples);
#'     "Gini" (Gini index).
#' @inheritParams core
#' @return A data.frame of samples x global indicators; except when split=FALSE, a vector of indices is returned.
#' @details This function returns the indices with the default choices for detection, prevalence and other parameters for simplicity and standardization. See the individual functions for more options. This function extends the functionality of the phyloseq::estimate_richness function.
#' @examples 
#'   data(dietswap)
#'   d <- global(dietswap)
#' @export
#' @seealso rarity, core_abundance, top_abundance, low_abundance, dominance, gini, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
global <- function(x, split = TRUE, measures = NULL) {

  res <- NULL
  
  selected.vegan <- c("Shannon", "InvSimpson")
  nonveg <- c("Richness", "Evenness", "Dominance", "Top_Abundance", "Gini", "Low_Abundance", "Core_Abundance")
  
  if (is.null(measures)) {
    measures <- unique(c(selected.vegan, nonveg))
  }

  # Remove non-vegan measures
  measures.veg <- setdiff(measures, nonveg)

  if (length(measures.veg) > 0) {
    res <- estimate_richness(x, split = split, measures = measures.veg)
  } 

  if (("Richness" %in% measures) || is.null(measures)) {

    ri <- richness(x, detection = 0, split)
    
    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Richness = ri)
    } else {
      res$Richness <- ri
    }

  }

  if (("Evenness" %in% measures) || is.null(measures)) {

    # Shannon Diversity
    if (!"Shannon" %in% names(res)) {
      d <- estimate_richness(x, split = split, measures = "Shannon")$Shannon
    } else {
      d <- res$Shannon
    }

    # Richness
    # Calculate here making sure detection = 0 as it is also for Shannon
    r <- richness(x, detection = 0, split)

    # normalize by Log richness to get Pielou's evenness
    e <- d/log(r)

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Evenness = e)
    } else {
      res$Evenness <- e
    }    
  }


  if (("Dominance" %in% measures) || is.null(measures)) {

    do <- unname(dominance(x, split = split))

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

  if (("Top_Abundance" %in% measures) || is.null(measures)) {

    do <- unname(top_abundance(x, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Top_Abundance = do)
    } else {
      res$Top_Abundance <- do
    }

  }

  if (("Low_Abundance" %in% measures) || is.null(measures)) {

    th <- quantile(as.vector(abundances(x)), 1)
    do <- unname(low_abundance(x, detection = 0.2/100, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Low_Abundance = do)
    } else {
      res$Low_Abundance <- do
    }

  }


  if (("Core_Abundance" %in% measures) || is.null(measures)) {

    do <- unname(core_abundance(x, detection = 0.2/100, prevalence = 50/100, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Core_Abundance = do)
    } else {
      res$Core_Abundance <- do
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





richness <- function (x, detection, split) {

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

    ri

}