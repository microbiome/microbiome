#' @title bimodality
#' @description Wrapper to calculate bimodality score based on various methods 
#'
#' @param x A vector, matrix, or a phyloseq object
#' @param method bimodality quantification method ('potential_bootstrap'
#' 	  or one of the methods in coefficient_of_bimodality)
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
#'     \item{Sarle.finite.sample}{Coefficient of bimodality for finite sample. See SAS 2012.}
#'     \item{Sarle.asymptotic}{Coefficient of bimodality, used and described in Shade et al. (2014) and Ellison AM (1987).}
#'     \item{potential_bootstrap}{Repeats potential analysis (Livina et al. 2010) multiple times with bootstrap sampling for each row of the input data (as in Lahti et al. 2014) and returns the bootstrap score.}
#'     \item{dip}{DIP test from the package \pkg{diptest}. See help(dip.test). The DIP test quantifies unimodality and varies between [0,1]. To quantify multimodality, this function returns the inverse score: 1-x. }
#'   }
#'
#' @seealso coefficient_of_bimodality
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
#' SAS Institute Inc. (2012). SAS/STAT 12.1 user's guide. Cary, NC.
#'
#' @importFrom diptest dip.test
#' @export
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples bimodality(c(rnorm(100, mean = 0), rnorm(100, mean = 5)))
#'
#' @keywords utilities
bimodality <- function (x, method = "potential_bootstrap", detection.threshold = 1, bw.adjust = 1, bs.iterations = 100, detection.limit = 1, verbose = TRUE) {

  if (is.vector(x)) {

    if (method %in% c("Sarle.finite.sample", "Sarle.asymptotic")) {

      s <- coefficient_of_bimodality(x, type = method)

    } else if (method == "potential_bootstrap") {

      if (length(unique(x)) == 1) {

        s <- 0

      } else {

        s <- multimodality_score(x, detection.threshold, 
      	   		       bw.adjust, bs.iterations, 
     			       detection.limit, verbose)$score
      }
    } else if (method == "dip") {

      # Pick OTU log10 data
      score <- dip.test(x, simulate.p.value = TRUE, B = 200)$statistic
      #data.frame(t(sapply(dip, function (x) {c(x$statistic, x$p.value)})))
      #colnames(dip2) <- c("score", "p.value")
      #dip2$tax <- names(dip)

      # Dip quantifies unimodality. Values range between 0 to 1. 
      # Values less than 0.05 indicate significant deviation from unimodality. 
      # To score multimodality, use the inverse:
      s <- 1 - score
    }

  } else if (is.matrix(x)) {

    s <- apply(x, 1, function (xi) {bimodality(xi, method = method, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit, verbose = verbose)})

  } else if (class(x) == "phyloseq") {

    # Pick the data from phyloseq object
    x <- log10(otu_table(x)@.Data)  
    s <- bimodality(x, method = method, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit, verbose = verbose)

  }
  
  s 

}


  
