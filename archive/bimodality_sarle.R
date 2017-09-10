#' @title Sarle's Bimodality Coefficient
#' @description Sarle's bimodality coefficient.
#' @param x Data vector for which bimodality will be quantified
#' @param bs.iter Bootstrap iterations
#' @param na.rm Remove NAs
#' @param type Score type ("Sarle.finite.sample" or "Sarle.asymptotic")
#' @return Bimodality score
#' @export
#' @examples
#'   b <- bimodality_sarle(rnorm(100), type = "Sarle.finite.sample")
#' @details The coefficient lies in (0, 1).
#' 
#'	    The 'Sarle.asymptotic' version is defined as
#'          \deqn{b = (g^2 + 1) / k}.
#'          This is coefficient of bimodality from Ellison AM Am. J. Bot. 1987, 
#'          for microbiome analysis it has been used for instance in
#'          Shade et al. 2014.
#'
#'          The formula for 'Sarle.finite.sample' (SAS 2012):
#'
#' 	    \deqn{b = \frac{g^2 + 1}{k + (3(n-1)^2)/((n-2)(n-3))}}
#'          where n is sample size and 
#' 
#'          In both formulas, \eqn{g} is sample skewness and \eqn{k} is the kth
#'          standardized moment (also called the sample kurtosis, or
#'          excess kurtosis).
#'
#' @references
#'   \itemize{
#'     \item{}{Shade et al. mBio 5(4):e01371-14, 2014.}
#'     \item{}{Ellison AM (1987) Am J Botany 74(8):1280-1288.}
#'     \item{}{SAS Institute Inc. (2012). SAS/STAT 12.1 user's guide. Cary, NC.}
#'     \item{}{To cite the microbiome R package, see citation('microbiome')}
#'  }
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso Check the dip.test from the \pkg{DIP} package for a
#' classical test of multimodality.
#' @keywords utilities
bimodality_sarle <- function(x, bs.iter = 1, na.rm = TRUE, type = "Sarle.finite.sample") {

    g <- skewness(x, na.rm)
    k <- kurtosis(x, na.rm)      

    if (type == "Sarle.asymptotic") {

      s <- (1 + g^2)/(k + 3)

    } else if (type == "Sarle.finite.sample") {
    
      n <- length(x)
      s <- (g^2 + 1) / (k + (3*(n-1)^2)/((n-2)*(n-3)))
      
    }

    if (bs.iter > 1) {
      s <- c()
      for (i in 1:bs.iter) {
        xbs <- sample(x, replace = TRUE)
	s[[i]] <- bimodality_sarle(xbs, type = type)
      }
      s <- mean(s)
    }

    s
    
}
