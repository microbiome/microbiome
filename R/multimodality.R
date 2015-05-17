#' multimodality
#'
#' Wrapper to calculate multimodality score based on various methods 
#'
#' @param x A vector, matrix, or a phyloseq object
#' @param method multimodality quantification method ('potential_bootstrap'
#' 	  or 'coefficient_of_bimodality')
#' @param detection.threshold Mode detection threshold
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iterations Bootstrap iterations
#' @param detection.limit minimum accepted density for a maximum; 
#'        as a multiple of kernel height
#' @param verbose Verbose
#' 
#' @return A list with following elements: 
#'   \itemize{
#'     \item{score}{Fraction of bootstrap samples where multiple modes are observed}
#'     \item{nmodes}{The most frequently observed number of modes in bootrstrap sampling results}
#'     \item{results}{Full results of potential_analysis_bootstrap for each row of the input matrix.}
#'   }
#'
#' @details
#'   \itemize{
#'     \item{coefficient_of_bimodality}{Coefficient of bimodality, used and described in Shade et al. (2014) and Ellison AM (1987).}
#'     \item{potential_bootstrap}{Repeats potential analysis (Livina et al. 2010) multiple times with bootstrap sampling for each row of the input data (as in Lahti et al. 2014) and returns the bootstrap score.}
#'   }
#'
#' @references
#' Livina et al. (2010). Potential analysis 
#' reveals changing number of climate states during the last 60
#' kyr. \emph{Climate of the Past}, 6, 77-82.
#'
#' Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.
#'
#' Shade et al. mBio 5(4):e01371-14, 2014.
#' 
#' AM Ellison, Am. J. Bot 74:1280-8, 1987.
#'
#' @export
#'
#' @import earlywarnings
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples multimodality(c(rnorm(100, mean = 0), rnorm(100, mean = 5)))
#'
#' @keywords utilities

multimodality <- function (x, method = "potential_bootstrap", detection.threshold = 1, bw.adjust = 1, bs.iterations = 100, detection.limit = 1, verbose = TRUE) {

  if (is.vector(x)) {

    if (method == "coefficient_of_bimodality") {

      s <- coefficient_of_bimodality(x)

    } else if (method == "potential_bootstrap") {


      if (length(unique(x)) == 1) {
        s <- 0
      } else {

        s <- multimodality_score(x, detection.threshold, 
      	   		       bw.adjust, bs.iterations, 
     			       detection.limit, verbose)$score
      }
    }

  } else if (is.matrix(x)) {

    s <- apply(x, 1, function (xi) {multimodality(xi, method = method, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit, verbose = verbose)})

  } else if (class(x) == "phyloseq") {

    # Pick the data from phyloseq object
    x <- log10(otu_table(x)@.Data)  
    s <- multimodality(x, method = method, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit, verbose = verbose)

  }
  
  s 

}


  
