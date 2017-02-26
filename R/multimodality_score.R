#' @title Multimodality Score
#' @description Multimodality score based on bootstrapped potential analysis.
#' @param x A vector, or data matrix (variables x samples)
#' @param detection Mode detection
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iterations Bootstrap iterations
#' @param detection.limit minimum accepted density for a maximum;
#'        as a multiple of kernel height
#' @param verbose Verbose
#' @return A list with following elements: 
#'   \itemize{
#'     \item{score}{Fraction of bootstrap samples with multiple observed modes}
#'     \item{nmodes}{The most frequently observed number of modes in bootstrap}
#'     \item{results}{Full results of potential_analysis_bootstrap for each
#'                      row of the input matrix.}
#' }
#' @details Repeats potential analysis (Livina et al. 2010) multiple times
#'          with bootstrap sampling for each row of the input data
#'          (as in Lahti et al. 2014) and returns the specified results.
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   data(peerj32)
#'   s <- multimodality_score(
#'       	    t(peerj32$microbes[, c("Akkermansia", "Dialister")]))
#' @references
#'  \itemize{
#'   \item{}{Livina et al. (2010). Potential analysis reveals changing number
#'             of climate states during the last 60 kyr.
#'	     \emph{Climate of the Past}, 6, 77-82.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#'     	     ecosystem. \emph{Nature Communications} 5:4344.}
#'  }
#' @keywords utilities
multimodality_score <- function (x, detection = 1, bw.adjust = 1, bs.iterations = 100, detection.limit = 1, verbose = TRUE) {

  if (is.vector(x)) {
  
    # Add small noise to enable robust density estimation
    # (identical values may cause failure)
    x <- x + rnorm(length(x), sd = sd(x)/100)
    m <- potential_analysis_bootstrap(x,
      	   detection = detection, bw.adjust = bw.adjust,
	   bs.iterations = bs.iterations, detection.limit = detection.limit)
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
			       bs.iterations, detection.limit, verbose)
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




