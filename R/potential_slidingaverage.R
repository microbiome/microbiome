#' @title Moving Average Potential
#' @description This function reconstructs a potential derived from data
#'              along a gradient of a given parameter.
#' @param X a vector of the X observations of the state variable of interest
#' @param param parameter values corresponding to the observations in X 
#' @param bw Bandwidth for smoothing kernels. Automatically determined
#'           by default.
#' @param bw.adjust Bandwidth adjustment constant
#' @param detection Threshold for local optima to be discarded.
#' @param std Standard deviation.
#' @param grid.size number of evaluation points; number of steps between
#'           min and max potential; also used as kernel window size
#' @param plot.cutoff cuttoff for potential minima and maxima in visualization
#' @param plot.contours Plot contours on the landscape visualization
#' @param binwidth binwidth for contour plot
#' @param bins bins for contour plot. Overrides binwidth if given
#' @return A list with the following elements:
#'   \itemize{
#'     \item{pars}{values of the covariate parameter as matrix}
#'     \item{xis}{values of the x as matrix}
#'     \item{pots}{smoothed potentials}
#'     \item{mins}{minima in the densities (-potentials; neglecting local optima)}
#'     \item{maxs}{maxima in densities (-potentials; neglecting local optima)}
#'     \item{plot}{an object that displays the potential estimated in 2D}
#'   }
#' @export
#' @references
#'  \itemize{
#'   \item{}{Hirota, M., Holmgren, M., van Nes, E.H. & Scheffer, M. (2011).
#'           Global resilience of tropical forest and savanna to critical
#'           transitions. \emph{Science}, 334, 232-235.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#'           ecosystem. \emph{Nature Communications} 5:4344.}
#'  }
#' @author Leo Lahti, adapted from original Matlab code by Egbert van Nes.
#' @seealso \code{potential_univariate}
#' @examples
#'   X <- c(rnorm(1000, mean = 0),
#'            rnorm(1000, mean = -2),
#'            rnorm(1000, mean = 2));
#'	      param = seq(0,5,length=3000); 
#'	    res <- potential_slidingaverage(X, param)
#' @keywords utils
potential_slidingaverage <- function(X, param = NULL, bw = "nrd",
                                     bw.adjust = 1, detection = 0.1,
				     std = 1, grid.size = 50,
				     plot.cutoff = 0.5,
				     plot.contours = TRUE, binwidth = 0.2, 
    				     bins = NULL) {
    
    if (is.null(param)) {
        param <- seq(1, length(X), 1)
    }
    
    nas <- is.na(param) | is.na(X)
    if (sum(nas) > 0) {
        warning("The data contains NAs, removing the associated samples from 
                 X and param input arguments.")
        X <- X[!nas]
        param <- param[!nas]
    }
    
    minparam <- min(param)
    maxparam <- max(param)
    
    # Determine step size
    sdwindow <- step <- (maxparam - minparam)/grid.size
    
    # Place evaluation points evenly across data range
    xi <- seq(min(X), max(X), length = grid.size)
    
    # Initialize
    xis <- pars <- pots <- matrix(0, nrow = grid.size, ncol = length(xi))
    maxs <- mins <- matrix(0, nrow = grid.size, ncol = length(xi))
    
    for (i in 1:grid.size) {
        
        # Increase the parameter at each step
        par <- minparam + (i - 0.5) * step
        
        # Check which elements in evaluation range (param) are within
	# 2*sd of par
        weights <- exp(-0.5 * (abs(par - param)/sdwindow)^2)
        
        # LL: Normalization was added in the R implementation 16.5.2012
        weights <- weights/sum(weights)
        
        # Calculate the potential
        tmp <- potential_univariate(x = X, std = std, bw = bw,
	                            bw.adjust = bw.adjust, 
            weights = weights, grid.size = grid.size)
        
        # Store variables
        pots[i, ] <- tmp$pot
        xis[i, ] <- tmp$grid.points
        pars[i, ] <- par + rep(0, length(tmp$grid.points))
        mins[i, tmp$min.inds] <- 1
        maxs[i, tmp$max.inds] <- 1
        
    }
    
    res <- list(pars = pars, xis = xis, pots = pots,
    	        mins = mins, maxs = maxs, std = std)
    
    p <- plot_potential(res, cutoff = plot.cutoff,
                        plot.contours = plot.contours,
			binwidth = binwidth, bins = bins)
			
    p <- p + xlab("parameter/time") + ylab("state variable")
        
    list(res = res, plot = p)
    
}

