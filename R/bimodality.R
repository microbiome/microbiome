#' @title coefficient_of_bimodality
#' @description Coefficient of bimodality (Sarle's bimodality coefficient b)
#' @param x Data vector for which bimodality will be quantified
#' @param bs.iter Bootstrap iterations
#' @param na.rm Remove NAs
#' @param type Bimodality score type ("Sarle.finite.sample" or "Sarle.asymptotic")
#' @return Bimodality score
#' @examples # coefficient_of_bimodality(rnorm(100), type = "Sarle.finite.sample")
#' @export
#' @importFrom moments skewness
#' @importFrom moments kurtosis
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
#' @keywords utilities
#' 
coefficient_of_bimodality <- function(x, bs.iter = 1, na.rm = TRUE, type = "Sarle.finite.sample") {

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
	s[[i]] <- coefficient_of_bimodality(xbs, type = type)
      }
      s <- mean(s)
    }

    s
    
}
