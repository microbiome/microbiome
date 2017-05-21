#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @param index If the index is given, it will override the other parameters.
#'     See the details below for description and references of the standard
#'     rarity indices. 
#' @inheritParams global
#' @return A vector of rarity indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rarity(dietswap, index = 'low_abundance')
#'   # d <- rarity(dietswap, index = 'all')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso global, log_modulo_skewness, noncore_abundance, low_abundance
#' @details
#'   The rarity index characterizes the concentration of species at
#'   low abundance.
#'
#'   The following rarity indices are provided:
#'   \itemize{
#'     \item{log_modulo_skewness}{Quantifies the concentration of the least
#'       abundant species by the log-modulo skewness of the arithmetic
#'       abundance classes (see Magurran & McGill 2011). These are typically
#'       right-skewed; to avoid taking log of occasional negative skews,
#'       we follow Locey & Lennon (2016) and use the log-modulo
#'       transformation that adds a value of one to each measure of skewness
#'       to allow logarithmization. The values q=0.5 and n = 50 are used here.}
#'     \item{low_abundance}{Relative proportion of the least abundant species,
#'        below the detection level of 0.2\%. The least abundant species are
#'        determined separately for each sample regardless of their prevalence.}
#'     \item{noncore_abundance}{Relative proportion of the non-core species,
#'        exceed the detection level of 0.2\% at 50\% prevalence at most.
#'        This is complement of the core with the same thresholds.}
#'     \item{rare_abundance}{Relative proportion of the rare 
#'        taxa in [0,1] - the rare taxa are detected with less than 20\%
#'        prevalence, regardless of abundance.}
#'   }
#' 
#' @references
#'   Kenneth J. Locey and Jay T. Lennon.
#'   Scaling laws predict global microbial diversity.
#'   PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#'   Magurran AE, McGill BJ, eds (2011)
#'   Biological Diversity: Frontiers in Measurement and Assessment
#'   (Oxford Univ Press, Oxford), Vol 12
#'
#' @keywords utilities
rarity <- function(x, index = "all") {
    
    # Only include accepted indices
    accepted <- c("log_modulo_skewness", "low_abundance",
                    "noncore_abundance", "rare_abundance")
    
    # Return all indices
    if ("all" %in% index) {
        index <- accepted
    }
    
    index <- intersect(index, accepted)
    if (length(index) == 0) {
        return(NULL)
    }
    
    
    if (length(index) > 1) {
        tab <- NULL
        for (idx in index) {
            tab <- cbind(tab, rarity(x, index = idx))
        }
        colnames(tab) <- index
        return(as.data.frame(tab))
    }
    
    # Pick data
    otu <- abundances(x)
    otu.relative <- abundances(transform(x, "compositional"))
    
    if (index == "log_modulo_skewness") {
    
        r <- log_modulo_skewness(otu, q = 0.5, n = 50)
	
    } else if (index == "low_abundance") {
    
        r <- apply(otu.relative, 2,
                function(x) low_abundance(x, detection = 0.2/100))
    } else if (index == "noncore_abundance") {
    
        r <- noncore_abundance(x, detection = 0.2/100,
                                prevalence = 20/100)
		
    } else if (index == "rare_abundance") {
    
        s <- rare_members(otu.relative,
                detection = 100/100, # All abundances accepted
                prevalence = 20/100) # Less than this prevalence required
        r <- colSums(otu.relative[s, ])
	
    }
    
    names(r) <- colnames(otu)
    
    r
    
}



