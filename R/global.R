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
#'     "richness" (number of unique taxa that give non-zero signal); 
#'     "evenness" (Pielou's index);
#'     "coverage" (Number of species needed to cover 50% of the ecosystem);
#'     "DBP" (Berger-Parker dominance ie. relative proportion of the single most dominant species in [0,1]); see ?dominance for details and references
#'     "DMN" (McNaughton dominance ie. relative proportion of the two most dominant species in [0,1]); see ?dominance for details and references
#'     "rarity" Quantifies the concentration of the least abundant species by the log-modulo skewness of the arithmetic abundance classes. See ?rarity for details and references.
#'     "low_abundance" (Relative proportion of the least abundant species, below the detection level of 0.2%); 
#'     "rare_abundance" (Relative proportion of the rare (non-core) species in [0,1]) - this complement (1-x) of the core_abundance
#'     "core_abundance" (Relative proportion of the core species that exceed detection level 0.2% in over 50% of the samples);
#'     "gini" (Gini index; calculated with the function inequality).
#' @inheritParams core
#' @return A data.frame of samples x global indicators; except when split=FALSE, a vector of indices is returned.
#' @details This function returns the indices with the default choices for detection, prevalence and other parameters for simplicity and standardization. See the individual functions for more options. This function extends the functionality of the phyloseq::estimate_richness function.
#' @examples 
#'   data(dietswap)
#'   d <- global(dietswap)
#' @export
#' @seealso core_abundance, dominance, rarity, coverage, gini, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
global <- function(x, split = TRUE, measures = NULL) {

  res <- NULL

  selected.vegan <- c("Shannon", "InvSimpson", "Fisher")
  
  nonveg <- c("richness", "evenness", "coverage", "DBP", "DMN", "gini", "rarity", "low_abundance", "core_abundance")
  
  if (is.null(measures)) {
    measures <- unique(c(selected.vegan, nonveg))
  }

  # Allow inputs also in full lowercase but
  # capitalize those for phyloseq::estimate_richness
  measures <- gsub("shannon", "Shannon", measures)
  measures <- gsub("invsimpson", "InvSimpson", measures)
  measures <- gsub("chao1", "Chao1", measures)
  measures <- gsub("ace", "ACE", measures)      

  # Remove non-vegan measures
  measures.veg <- setdiff(measures, nonveg)

  if (length(measures.veg) > 0) {
    res <- estimate_richness(x, split = split, measures = measures.veg)
  } 

  if (("richness" %in% measures) || is.null(measures)) {

    ri <- richness(x, detection = 0, split)
    
    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(richness = ri)
    } else {
      res$richness <- ri
    }

  }

  if (("evenness" %in% measures) || is.null(measures)) {

    # Shannon Diversity
    if (!"Shannon" %in% names(res)) {
      d <- estimate_richness(x, split = split, measures = "Shannon")$Shannon
    } else {
      d <- res$Shannon
    }

    # richness
    # Calculate here making sure detection = 0 as it is also for Shannon
    r <- richness(x, detection = 0, split)

    # normalize by Log richness to get Pielou's evenness
    e <- d/log(r)

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(evenness = e)
    } else {
      res$evenness <- e
    }    
  }


  if (("coverage" %in% measures) || is.null(measures)) {

    do <- unname(coverage(x, split = split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(coverage = do)
    } else {
      res$coverage <- do
    }

  }


  if (("gini" %in% measures) || is.null(measures)) {

    do <- unname(inequality(x, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(gini = do)
    } else {
      res$gini <- do
    }

  }

  if (any(c("dominance", "DBP") %in% measures) || is.null(measures)) {

    do <- unname(dominance(x, index = "DBP", split = TRUE))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(DBP = do)
    } else {
      res$DBP <- do
    }
  }

  if ((c("DMN") %in% measures) || is.null(measures)) {

    do <- unname(dominance(x, index = "DMN", split = TRUE))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(DMN = do)
    } else {
      res$DMN <- do
    }
  }

  if (("rarity" %in% measures) || is.null(measures)) {

    if (!split) {
      stop("split not defined for the rarity index")
    }

    do <- unname(rarity(x))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(rarity = do)
    } else {
      res$rarity <- do
    }
  }

  if (("core_abundance" %in% measures) || is.null(measures)) {

    do <- unname(core_abundance(x, detection = 0.2/100, prevalence = 50/100, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(core_abundance = do)
    } else {
      res$core_abundance <- do
    }
  }

  if (("rare_abundance" %in% measures) || is.null(measures)) {

    do <- unname(rare_abundance(x, detection = 0.2/100, prevalence = 50/100, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(rare_abundance = do)
    } else {
      res$rare_abundance <- do
    }

  }

  if (("low_abundance" %in% measures) || is.null(measures)) {

    do <- unname(low_abundance(x, detection = 0.2/100, split))

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(low_abundance = do)
    } else {
      res$low_abundance <- do
    }
  }


  # For total diversity, just return a vector of values
  # (as the table would have only 1 row; and giving a different
  #  output for this special case can potentially help avoid confusion)
  if (!split) {
    res <- unlist(res)
  }

  # All indicator names in lowercase for clarity and consistency
  if (length(res) > 1) { 
    colnames(res) <- tolower(colnames(res))
  } else {
    names(res) <- tolower(names(res))
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