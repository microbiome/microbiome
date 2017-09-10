#' @title Phyloseq Multimodality Test
#' @description Calculate multimodality score based on bootstrapped
#' potential analysis.
#' @param x Phyloseq object
#' @param detection Mode detection
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iterations Bootstrap iterations
#' @param detection.limit minimum accepted density for a maximum;
#'        as a multiple of kernel height
#' @param verbose Verbose
#' @return A list with following elements
#'   \item{score}{Fraction of bootstrap samples where multiple modes
#'                are observed}
#'   \item{nmodes}{The most frequently observed number of modes in bootstrap
#'                 sampling results}
#'   \item{results}{Full results of potential_analysis_bootstrap for each
#'                  row of the input matrix}
#' @details Repeats potential analysis (Livina et al. 2010) multiple times
#' with bootstrap sampling for each row of the input data
#'  (as in Lahti et al. 2014) and returns the specified results.
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   data(peerj32)
#'   s <- multimodality_score(
#'            t(peerj32$microbes[, c("Akkermansia", "Dialister")]))
#' @references
#'  \itemize{
#'   \item{}{Livina et al. (2010). Potential analysis reveals changing number
#'           of climate states during the last 60 kyr.
#'           \emph{Climate of the Past}, 6, 77-82.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#'           ecosystem. \emph{Nature Communications} 5:4344.}
#'  }
#' @keywords utilities
multimodality_phyloseq <- function (x, detection = 1, bw.adjust = 1,
                                    bs.iterations = 100, detection.limit = 1,
				    verbose = TRUE) {
  x <- abundances(transform_phyloseq(x, "log10"))
  msc <- multimodality_score(x, detection,
      	   bw.adjust, bs.iterations, detection.limit, verbose)
  msc
  
}

