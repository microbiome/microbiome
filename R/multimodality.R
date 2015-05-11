#' multimodality
#'
#' Wrapper to calculate multimodality score based on various methods 
#'
#' @param x A vector, matrix, or a phyloseq object
#' @param method multimodality quantification method ("potential.bootstrap" or "coefficient.of.bimodality")
#' @param detection.threshold Mode detection threshold
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iterations Bootstrap iterations
#' @param detection.limit minimum accepted density for a maximum; as a multiple of kernel height
#' @param verbose Verbose
#' 
#' @return A list with following elements: 
#' 	  \itemize{
#'		\item{score}{Fraction of bootstrap samples where multiple modes are observed}
#'	     	\item{nmodes}{The most frequently observed number of modes in bootrstrap sampling results}
#'		\item{results}{Full results of potential_analysis_bootstrap for each row of the input matrix.}
#'	   }
#'
#' @details
#' 	  \itemize{
#'	    \item{coefficient.of.bimodality}{Coefficient of bimodality, used and described in Shade et al.
#'	    				     (2014) and Ellison AM (1987).}
#'  	    \item{potential.bootstrap}{Repeats potential analysis (Livina et al. 2010) multiple times with bootstrap
#'	    			sampling for each row of the input data (as in Lahti et al. 2014) and returns the
#'				bootstrap score.}
#'       }
#'
#' @export
#'
#' @import earlywarnings
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   data(peerj32)
#'   multimodality(peerj32$microbes[, "Akkermansia"])
#'
#' @references
#' Livina et al. (2010). Potential analysis 
#' reveals changing number of climate states during the last 60
#' kyr. \emph{Climate of the Past}, 6, 77-82.
#'
#' Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.
#'
#' Shade et al. (2014). Conditionally Rare Taxa Disproportionately Contribute to Temporal Changes in Microbial Diversity.
#' mBio 5(4):e01371-14. doi: 10.1128/mBio.01371-14
#'
#' Ellison AM (1987). Effect of seed dimorphism on the density-dependent dynamics of experimental populations of
#' Atriplex triangularis (Chenopodiaceae). Am. J. Bot. 74:1280–1288. doi:10.2307/2444163.
#'
#' @keywords utilities

multimodality <- function (x, method = "potential.bootstrap", detection.threshold = 1, bw.adjust = 1, bs.iterations = 100, detection.limit = 1, verbose = TRUE) {

  if (is.vector(x)) {
    if (method == "coefficient.of.bimodality") {
      s <- coefficient.of.bimodality(x)
    } else if (method == "potential.bootstrap") {
      s <- multimodality_score(x, detection.threshold, bw.adjust, bs.iterations, detection.limit, verbose)$score
    }
  } else if (is.matrix(x)) {
    s <- apply(x, 1, function (x) {multimodality(x, method = method, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit, verbose = verbose)})
  } else if (class(x) == "phyloseq") {
    # Pick the data from phyloseq object
    x <- log10(otu_table(x)@.Data)  
    s <- apply(x, 1, function (x) {multimodality(x, method = method, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit, verbose = verbose)}) 
  }
  
  s 

}


  
