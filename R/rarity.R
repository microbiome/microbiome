#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @inheritParams core
#' @param q Arithmetic abundance classes are evenly cut up to to this quantile of the data. The assumption is that abundances higher than this are not common, and they are classified in their own group.
#' @param n The number of arithmetic abundance classes from zero to the quantile cutoff indicated by q.
#' @return A vector of rarity indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rarity(dietswap)
#' @details The rarity index characterizes the concentration of species at low abundance. Here, we use the skewness of the frequency distribution of arithmetic abundance classes (see Magurran & McGill 2011). These are typically right-skewed; to avoid taking log of occasional negative skews, we follow Locey & Lennon (2016) and use the log-modulo transformation that adds a value of one to each measure of skewness to allow logarithmization. 
#' @references
#'   Kenneth J. Locey and Jay T. Lennon. Scaling laws predict global microbial diversity. PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#'   Magurran AE, McGill BJ, eds (2011) Biological Diversity: Frontiers in Measurement and Assessment (Oxford Univ Press, Oxford), Vol 12
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, low_abundance, global
rarity <- function (x, q = .5, n = 50) {

  a <- abundances(x)

  # Determine the quantile point. 
  th1 <- quantile(max(a), q)

  # Tabulate the arithmetic abundance classes
  # Use the same classes for all samples for consistency
  tab <- table(cut(as.vector(a), c(seq(0, th1, length = n), Inf)))

  # Check skewness of the abundance classes for each sample
  r <- apply(a, 2, function (x) { skew(tab) })

  # Return log-modulo
  log(1 + r)

}


# Inspired by moments::skewness but rewritten.
# Internal.
skew <- function (x) {
  x <- x[!is.na(x)]
  n <- length(x)
  (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
}