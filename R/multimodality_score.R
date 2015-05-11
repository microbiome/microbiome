#' multimodality_phyloseq
#'
#' Calculate multimodality score based on bootstrapped potential analysis
#'
#' @param x Phyloseq object
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
#' @details This function repeats potential analysis (Livina et al. 2010) multiple 
#' 	    times with bootstrap sampling for each row of the input data 
#'	    (as in Lahti et al. 2014) and returns the specified results.
#'
#' @export
#'
#' @import earlywarnings
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   data(peerj32)
#'   multimodality_score(t(peerj32$microbes[, c("Akkermansia", "Dialister")]))
#'
#' @references 
#' Livina et al. (2010). Potential analysis 
#' reveals changing number of climate states during the last 60
#' kyr. \emph{Climate of the Past}, 6, 77-82.
#'
#' Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.
#'
#' @keywords utilities

multimodality_phyloseq <- function (x, detection.threshold = 1, bw.adjust = 1, bs.iterations = 100, detection.limit = 1, verbose = TRUE) {

  x <- log10(otu_table(x)@.Data)
  msc <- multimodality_score(x, detection.threshold, bw.adjust, bs.iterations, detection.limit, verbose)
  msc
  
}

#' multimodality_score
#'
#' Calculate multimodality score based on bootstrapped potential analysis
#'
#' @param x A vector, or data matrix (variables x samples)
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
#' @details This function repeats potential analysis (Livina et al. 2010) multiple 
#' 	    times with bootstrap sampling for each row of the input data 
#'	    (as in Lahti et al. 2014) and returns the specified results.
#'
#' @export
#'
#' @import earlywarnings
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   data(peerj32)
#'   multimodality_score(t(peerj32$microbes[, c("Akkermansia", "Dialister")]))
#'
#' @references 
#' Livina et al. (2010). Potential analysis 
#' reveals changing number of climate states during the last 60
#' kyr. \emph{Climate of the Past}, 6, 77-82.
#'
#' Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.
#'
#' @keywords utilities
multimodality_score <- function (x, detection.threshold = 1, bw.adjust = 1, bs.iterations = 100, detection.limit = 1, verbose = TRUE) {

  if (is.vector(x)) {
  
    # Add small noise to enable robust density estimation (identical values may cause failure)
    x <- x + rnorm(length(x), sd = sd(x)/100)
    m <- potential_analysis_bootstrap(x, detection.threshold = detection.threshold, bw.adjust = bw.adjust, bs.iterations = bs.iterations, detection.limit = detection.limit)
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
      m <- multimodality_score(as.numeric(x[tax, ]), detection.threshold, bw.adjust, bs.iterations, detection.limit, verbose)
      nmodes[[tax]] <- m$modes 
      potential.results[[tax]] <- m
    }

    multimodality.score <- sapply(potential.results, function (x) { 1 - x$unimodality.support })
    ret <- list(score = multimodality.score, modes = nmodes, results = potential.results)

  }

  ret

}





#' potential_analysis_bootstrap
#'
#' Description:
#'
#' Bootstrap analysis of multimodality based on potential analysis of
#' Livina et al. (2010) as described in Lahti et al. (2014)
#' 
#' @param x Data vector
#' @param detection.threshold Mode detection threshold
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iterations Bootstrap iterations
#' @param detection.limit minimum accepted density for a maximum; as a multiple of kernel height
#'
#' @return List with following elements:
#' modes:  Number of modes for the input data vector (the most frequent number of modes from bootstrap)
#' minima: Average of potential minima across the bootstrap samples (for the most frequent number of modes)
#' maxima: Average of potential maxima across the bootstrap samples (for the most frequent number of modes)
#' unimodality.support Fraction of bootstrap samples exhibiting unimodality	   
#'
#' @import earlywarnings
#' @export
#' @references 
#' Livina et al. (2010). Potential analysis 
#' reveals changing number of climate states during the last 60
#' kyr. \emph{Climate of the Past}, 6, 77-82.
#'
#' Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.
#'
potential_analysis_bootstrap <- function (x, detection.threshold, bw.adjust = 1, bs.iterations = 100, detection.limit = 1) {

  nmodes <- c()
  minpoints <- list()
  maxpoints <- list()
  bws <- c()
  #s <- list()
  for (r in 1:bs.iterations) {
  
    # Bootstrap
    rs <- sample(length(x), replace = TRUE) 

    xbs <- na.omit(unname(x[rs]))
    #s[[r]] <- xbs

    a <- livpotential_ews(xbs, grid.size = floor(.2*length(x)), 
      	 		     detection.threshold = detection.threshold, 
			     bw.adjust = bw.adjust, 
			     detection.limit = detection.limit)

    nmodes[[r]] <- length(a$max.points)
    minpoints[[r]] <- a$min.points
    maxpoints[[r]] <- a$max.points
    bws[[r]] <- a$bw
  }

  # Most frequently observed number of modes
  top.modes <- as.numeric(names(which.max(table(nmodes))))
  min.points <- colMeans(do.call("rbind", minpoints[nmodes == top.modes]))
  max.points <- colMeans(do.call("rbind", maxpoints[nmodes == top.modes]))
  unimodality.support <- mean(nmodes <= 1)

  # Return the most frequent number of modes and
  # the corresponding tipping points
  # from the bootstrap analysis
  list(modes = top.modes, minima = min.points, maxima = max.points, unimodality.support = unimodality.support, bws = bws)
  
}

