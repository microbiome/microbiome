#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @param index If the index is given, it will override the other parameters.
#' See the details below for description and references of the standard
#' rarity indices. 
#' @inheritParams core
#' @return A vector of rarity indices
#' @export
#' @examples
#' data(dietswap)
#' d <- rarity(dietswap, index='low_abundance')
#' # d <- rarity(dietswap, index='all')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso alpha, log_modulo_skewness, rare_abundance, low_abundance
#' @details
#' The rarity index characterizes the concentration of species at
#' low abundance.
#'
#' The following rarity indices are provided:
#' \itemize{
#' \item{log_modulo_skewness }{Quantifies the concentration of the least
#' abundant species by the log-modulo skewness of the arithmetic
#' abundance classes (see Magurran & McGill 2011). These are typically
#' right-skewed; to avoid taking log of occasional negative skews,
#' we follow Locey & Lennon (2016) and use the log-modulo
#' transformation that adds a value of one to each measure of skewness
#' to allow logarithmization. The values q=0.5 and n=50 are used here.}
#' \item{low_abundance }{Relative proportion of the least abundant species,
#' below the detection level of 0.2\%. The least abundant species are
#' determined separately for each sample regardless of their prevalence.}
#' \item{rare_abundance }{Relative proportion of the non-core species,
#' exceed the given detection level (default 20%)
#' at the given prevalence (default 20%).
#' This is complement of the core with the same thresholds.}
#' }
#' 
#' @references
#' Kenneth J. Locey and Jay T. Lennon.
#' Scaling laws predict global microbial diversity.
#' PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#' Magurran AE, McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment
#' (Oxford Univ Press, Oxford), Vol 12
#'
#' @keywords utilities
rarity <- function(x, index = "all", detection = 0.2/100, prevalence = 20/100) {

    # Only include accepted indices
    index <- tolower(index)    
    accepted <- c("log_modulo_skewness", "low_abundance",
                "rare_abundance")
    accepted <- tolower(accepted)

    # Return all indices
    if (length(index) == 1 && index == "all") {
        index <- accepted
    }
    
    if (!is.null(index)) {
        index <- intersect(index, accepted)
    }
    
    if (!is.null(index) && length(index) == 0) {
        return(NULL)
    }

    tab <- rarity_help(x, index, detection, prevalence)

    if (is.vector(tab)) {
        tab <- as.matrix(tab, ncol=1)
        colnames(tab) <- index        
    }
    
    as.data.frame(tab) 

}


rarity_help <- function(x, index="all", detection, prevalence) {

    if ( length(index) > 1 ) {
    
        tab <- NULL
    
        for (idx in index) {
            tab <- cbind(tab,
        rarity_help(x, index = idx, detection, prevalence))
        }

        colnames(tab) <- index
        return(as.data.frame(tab))
    
    }

    # Pick data
    otu  <- abundances(x)
    otuc <- transform(x, "compositional")
    otu.relative <- abundances(otuc)

    if (index == "log_modulo_skewness") {
    
        r <- log_modulo_skewness(otu, q=0.5, n=50)
    
    } else if (index == "low_abundance") {

        r <- apply(otu.relative, 2,
                function(x) low_abundance(x, detection=detection))
        
    } else if (index == "rare_abundance") {
    
        r <- rare_abundance(x, detection=detection,
                            prevalence=prevalence)
            
    }
    
    names(r) <- colnames(otu)
    
    r
    
}



