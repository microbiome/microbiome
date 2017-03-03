#' @title Bimodality Analysis
#' @description A wrapper to calculate bimodality scores.
#' @param x A vector, matrix, or a phyloseq object
#' @param method bimodality quantification method ('potential_analysis'
#' 	  or one of the methods in the function \code{\link{bimodality_sarle}})
#' @param detection Mode detection
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iter Bootstrap iterations
#' @param min.density minimum accepted density for a maximum; as a multiple of kernel height
#' @param verbose Verbose
#' @return A list with following elements:
#'   \itemize{
#'     \item{score}{Fraction of bootstrap samples where multiple modes are observed}
#'     \item{nmodes}{The most frequently observed number of modes in bootrstrap sampling results}
#'     \item{results}{Full results of potential_analysis for each row of the input matrix.}
#'   }
#' @details
#'   \itemize{
#'     \item{Sarle.finite.sample}{Coefficient of bimodality for finite sample. See SAS 2012.}
#'     \item{Sarle.asymptotic}{Coefficient of bimodality, used and described in Shade et al. (2014) and Ellison AM (1987).}
#'     \item{potential_analysis}{Repeats potential analysis (Livina et al. 2010) multiple times with bootstrap sampling for each row of the input data (as in Lahti et al. 2014) and returns the bootstrap score.}
#'   }
#' @seealso A classical test of multimodality is provided by dip.test in the \pkg{DIP} package.
#' @references
#' \itemize{
#'   \item{}{Livina et al. (2010). Potential analysis 
#'         reveals changing number of climate states during the last 60
#' 	   kyr. \emph{Climate of the Past}, 6, 77-82.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#'         ecosystem. \emph{Nature Communications} 5:4344.}
#'   \item{}{Shade et al. mBio 5(4):e01371-14, 2014.}
#'   \item{}{AM Ellison, Am. J. Bot 74:1280-8, 1987.}
#'   \item{}{SAS Institute Inc. (2012). SAS/STAT 12.1 user's guide. Cary, NC.}
#' }
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#'   
#'   bimodality(c(rnorm(100, mean = 0), rnorm(100, mean = 5)))
#'  
#'  # The classical DIP test:
#'  # quantifies unimodality. Values range between 0 to 1. 
#'  # dip.test(x, simulate.p.value = TRUE, B = 200)$statistic
#'  # Values less than 0.05 indicate significant deviation from unimodality. 
#' @keywords utilities
bimodality <- function (x, method = "potential_analysis", detection = 1, bw.adjust = 1, bs.iter = 100, min.density = 1, verbose = TRUE) {

  if (is.vector(x)) {

    if (method %in% c("Sarle.finite.sample", "Sarle.asymptotic")) {

      s <- bimodality_sarle(x, type = method)

    } else if (method == "potential_analysis") {

      if (length(unique(x)) == 1) {

        s <- 0

      } else {

        # Shift the data. This does not affect mode detection but
	# avoids errors with nonnegatives.
        s <- multimodality_score(x, detection, 
      	   		       bw.adjust, bs.iter, 
     			       min.density, verbose)$score
      }
    }

  } else if (is.matrix(x)) {

    s <- apply(x, 1, function (xi) {bimodality(xi, method = method, detection = detection, bw.adjust = bw.adjust, bs.iter = bs.iter, min.density = min.density, verbose = verbose)})

  } else if (class(x) == "phyloseq") {

    # Pick the data from phyloseq object
    x <- abundances(x)
    s <- bimodality(x, method = method, detection = detection, bw.adjust = bw.adjust, bs.iter = bs.iter, min.density = min.density, verbose = verbose)

  }
  
  s 

}


  


#' @title Multimodality Score
#' @description Multimodality score based on bootstrapped potential analysis.
#' @param x A vector, or data matrix (variables x samples)
#' @param detection Mode detection
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iter Bootstrap iterations
#' @param min.density minimum accepted density for a maximum;
#'        as a multiple of kernel height
#' @param verbose Verbose
#' @return A list with following elements: 
#'   \itemize{
#'     \item{score}{Fraction of bootstrap samples with multiple observed modes}
#'     \item{nmodes}{The most frequently observed number of modes in bootstrap}
#'     \item{results}{Full results of potential_analysis for each
#'                      row of the input matrix.}
#' }
#' @details Repeats potential analysis (Livina et al. 2010) multiple times
#'          with bootstrap sampling for each row of the input data
#'          (as in Lahti et al. 2014) and returns the specified results.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#'   \dontrun{
#'     # Not exported
#'     data(peerj32)
#'     s <- multimodality_score(
#'       	    t(peerj32$microbes[, c("Akkermansia", "Dialister")]))
#'   }
#' @references
#'   \itemize{
#'     \item{}{Livina et al. (2010). Potential analysis reveals changing number
#'             of climate states during the last 60 kyr.
#'	       \emph{Climate of the Past}, 6, 77-82.}
#'     \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#'     	       ecosystem. \emph{Nature Communications} 5:4344.}
#'   }
#' @keywords utilities
multimodality_score <- function (x, detection = 1, bw.adjust = 1, bs.iter = 100, min.density = 1, verbose = TRUE) {

  if (is.vector(x)) {
  
    # Add small noise to enable robust density estimation
    # (identical values may cause failure)
    x <- x + rnorm(length(x), sd = sd(x)/100)
    m <- potential_analysis(x,
      	   bw.adjust = bw.adjust,
	   bs.iter = bs.iter,
	   min.density = min.density)
    ret <- list(score = 1 - m$unimodality.support, modes = m$modes, results = m)
    return(ret)

  } else {

    # Univariate potential analysis for all taxa with full data
    potential.results <- list()
    nmodes <- c()
    if (is.null(rownames(x))) {
      rownames(x) <- as.character(1:nrow(x))
    }

    for (tax in rownames(x)) {
      if (verbose) { message(tax) }
      m <- multimodality_score(as.numeric(x[tax, ]),
      	   		       detection, bw.adjust,
			       bs.iter, min.density, verbose)
      nmodes[[tax]] <- m$modes 
      potential.results[[tax]] <- m
    }

    multimodality.score <- sapply(potential.results,
    			function (x) { 1 - x$unimodality.support })
			
    ret <- list(score = multimodality.score,
    	        modes = nmodes,
		results = potential.results)

  }

  ret

}


