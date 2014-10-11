#' Coefficient of bimodality 
#'
#' @param x Data vector for which bimodality will be quantified
#'
#' @return Bimodality score
#'
#' @examples coefficient.of.bimodality(rnorm(100))
#' @export
#' @import e1071
#' @details Coefficient of bimodality used in Shade et al. mBio 5(4):e01371-14. 
#'          and picked from Ellison AM Am. J. Bot 74:1280-8, 1987 
#' @references 
#'   Shade et al. mBio 5(4):e01371-14, 2014.
#'   AM Ellison, Am. J. Bot 74:1280-8, 1987.
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
coefficient.of.bimodality <- function(x) {
    (1 + skewness(x)^2)/(kurtosis(x) + 3)
}
