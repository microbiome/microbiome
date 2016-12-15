#' @title Potential Analysis for Univariate Data
#' @description One-dimensional potential estimation for univariate timeseries.
#' @param x Univariate data (vector) for which the potentials shall be estimated
#' @param std Standard deviation of the noise (defaults to 1; this will set scaled potentials)
#' @param bw kernel bandwidth estimation method 
#' @param weights optional weights in ksdensity (used by
#' potential_slidingaverages).
#' @param grid.size Grid size for potential estimation.
#' @param detection.threshold maximum detection threshold as fraction
#' of density kernel height dnorm(0, sd = bandwidth)/N
#' @param bw.adjust The real bandwidth will be bw.adjust*bw; defaults to 1
#' @param density.smoothing Add a small constant density across the
#' whole observation range to regularize density estimation (and to
#' avoid zero probabilities within the observation range). This
#' parameter adds uniform density across the observation range, scaled
#' by density.smoothing.
#' @param detection.limit minimum accepted density for a maximum; as a
#' multiple of kernel height
#' @return \code{potential_univariate} returns a list with the
#' following elements:
#'   \itemize{
#'     \item{xi}{the grid of points on which the potential is estimated}
#'     \item{pot}{The estimated potential: -log(f)*std^2/2, where f is the density.}
#'     \item{density}{Density estimate corresponding to the potential.}
#'     \item{min.inds}{indices of the grid points at which the density has minimum values; (-potentials; neglecting local optima)}
#'     \item{max.inds}{indices the grid points at which the density has maximum values; (-potentials; neglecting local optima)}
#'     \item{bw}{bandwidth of kernel used}
#'     \item{min.points}{grid point values at which the density has minimum values; (-potentials; neglecting local optima)}
#'     \item{max.points}{grid point values at which the density has maximum values; (-potentials; neglecting local optima)}
#'  }
#' @export
#' @references
#'  \itemize{
#'   \item{}{Livina et al. (2010). Potential analysis reveals changing number of climate states during the last 60 kyr. \emph{Climate of the Past}, 6, 77-82.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal ecosystem. \emph{Nature Communications} 5:4344.}
#'  }
#' @author Based on Matlab code from Egbert van Nes modified by Leo Lahti. Extended from the initial version in the \pkg{earlywarnings} R package.
#' @seealso \code{\link{potential_slidingaverage}}
#' @examples \dontrun{res <- potential_univariate(x)}
#' @keywords early-warning
potential_univariate <- function(x, std = 1, bw = "nrd", weights = c(),
		     grid.size = NULL, 
    		     detection.threshold = 1, bw.adjust = 1,
		     density.smoothing = 0, detection.limit = 1) {
    
    if (is.null(grid.size)) {
        grid.size <- floor(0.2 * length(x))
    }
    
    # Density estimation
    tmp = try(de <- density(x, bw = bw, adjust = bw.adjust, 
       	  	  kernel = "gaussian", weights = weights, 
        	  window = kernel, n = grid.size, 
		  from = min(x), to = max(x), 
		  cut = 3, na.rm = FALSE))
    if (class(tmp) == "try-error") {
      # Just use default parameters if failing otherwise
      warning("Density estimation with custom parameters failed. 
      		       Using the defaults.")
      de <- density(x)
    }

    # Smooth the estimated density (f <- de$y) by adding a small
    # probability across the whole observation range (to avoid zero
    # probabilities for points in the observation range)
    f <- de$y + density.smoothing * 1/diff(range(de$x)) # *max(de$y)
    
    # Normalize the density such that it integrates to unity
    f <- f/sum(diff(de$x[1:2]) * f)
    
    # Final grid points and bandwidth
    grid.points <- de$x
    bw <- de$bw
    
    # Compute potential
    U <- -log(f) * std^2/2

    # Backtransform to density distribution
    # f <- exp(-2*U/std^2) 

    # Ignore very local optima

    # Note mins and maxs for density given # here (not for potential,
    # which has the opposite signs)
    ops <- find_optima(f, detection.threshold = detection.threshold, bw = bw, detection.limit = detection.limit)
    min.points <- grid.points[ops$min]
    max.points <- grid.points[ops$max]
    det.th <- ops$detection.threshold
    
    list(grid.points = grid.points, pot = U, density = f, min.inds = ops$min,
         max.inds = ops$max, bw = bw,
	 min.points = min.points, max.points = max.points,
	 detection.threshold = det.th)
    
}

