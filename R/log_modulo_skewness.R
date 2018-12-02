#' @title Log-Modulo Skewness Rarity Index
#' @description Calculates the community rarity index by log-modulo skewness.
#' @param x Abundance matrix (taxa x samples) with counts
#' @param q Arithmetic abundance classes are evenly cut up to to this quantile
#' of the data. The assumption is that abundances higher than this are not
#' common, and they are classified in their own group.
#' @param n The number of arithmetic abundance classes from zero to the
#' quantile cutoff indicated by q.
#' @return A vector of rarity indices
#' @examples
#' data(dietswap)
#' d <- log_modulo_skewness(dietswap)
#' @details The rarity index characterizes the concentration of species at low
#' abundance. Here, we use the skewness of the frequency distribution of
#' arithmetic abundance classes (see Magurran & McGill 2011).
#' These are typically right-skewed; to avoid taking log of occasional
#' negative skews, we follow Locey & Lennon (2016) and use the log-modulo
#' transformation that adds a value of one to each measure of skewness to
#' allow logarithmization. 
#' @references
#' Kenneth J. Locey and Jay T. Lennon.
#' Scaling laws predict global microbial diversity.
#' PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#' Magurran AE, McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment
#' (Oxford Univ Press, Oxford), Vol 12
#'
#' @export
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, low_abundance, alpha
log_modulo_skewness <- function(x, q=0.5, n=50) {

    # Get taxa x samples matrix
    a <- abundances(x)
    
    # Determine the quantile point.
    th1 <- quantile(max(a), q)
    
    # Tabulate the arithmetic abundance classes Use the same classes
    # for all samples for consistency    
    cutpoints <- c(seq(0, th1, length=n), Inf)
    
    # Check skewness of the abundance classes for each sample
    r <- apply(a, 2, function(x) {
        skew(table(cut(x, cutpoints)))
    })
    
    # Return log-modulo
    log(1 + r)
    
}


# Inspired by moments::skewness but rewritten.  Internal.
skew <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
}

