#' @title Bootstrapped Potential Analysis 
#' @description Analysis of multimodality based on bootstrapped potential
#'    analysis of Livina et al. (2010) as described in Lahti et al. (2014).
#' @param x Input data vector
#' @param detection.threshold Mode detection threshold
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iterations Bootstrap iterations
#' @param detection.limit minimum accepted density for a maximum; as a multiple of kernel height
#' @return List with following elements:
#' \itemize{
#'   \item{modes}{Number of modes for the input data vector (the most frequent number of modes from bootstrap)}
#'   \item{modes}{minima: Average of potential minima across the bootstrap samples (for the most frequent number of modes)}
#'   \item{modes}{maxima: Average of potential maxima across the bootstrap samples (for the most frequent number of modes)}
#'   \item{modes}{unimodality.support Fraction of bootstrap samples exhibiting unimodality}
#' }
#' @export
#' @references
#'  \itemize{
#'   \item{}{Livina et al. (2010). Potential analysis reveals changing number of climate states during the last 60 kyr. \emph{Climate of the Past}, 6, 77-82.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal ecosystem. \emph{Nature Communications} 5:4344.}
#'  }
potential_analysis_bootstrap <- function (x, detection.threshold, bw.adjust = 1, bs.iterations = 100, detection.limit = 1) {

  nmodes <- c()
  minpoints <- list()
  maxpoints <- list()
  bws <- c()

  for (r in 1:bs.iterations) {
  
    # Bootstrap
    rs <- sample(length(x), replace = TRUE) 

    xbs <- na.omit(unname(x[rs]))

    a <- potential_univariate(xbs, grid.size = floor(.2*length(x)), 
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
  list(modes = top.modes, minima = min.points, maxima = max.points,
       unimodality.support = unimodality.support, bws = bws)
  
}

