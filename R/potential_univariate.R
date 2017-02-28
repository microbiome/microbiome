#' @title Potential Analysis for Univariate Data
#' @description One-dimensional potential estimation for univariate timeseries.
#' @param x Univariate data (vector) for which the potentials shall be estimated
#' @param std Standard deviation of the noise (defaults to 1; this will set scaled potentials)
#' @param bw kernel bandwidth estimation method 
#' @param weights optional weights in ksdensity (used by
#' potential_slidingaverages).
#' @param grid.size Grid size for potential estimation.
#' @param detection maximum detection as fraction
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
    		     detection = 1, bw.adjust = 1,
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
    ops <- find_optima(f, detection = detection, bw = bw, detection.limit = detection.limit)
    min.points <- grid.points[ops$min]
    max.points <- grid.points[ops$max]
    det.th <- ops$detection
    
    list(grid.points = grid.points, pot = U, density = f, min.inds = ops$min,
         max.inds = ops$max, bw = bw,
	 min.points = min.points, max.points = max.points,
	 detection = det.th)
    
}





#' @title Find Optima
#' @description Detect optima, excluding local optima below
#'              detection. 
#' @param f density
#' @param detection detection for peaks
#' @param bw bandwidth
#' @param detection.limit Minimun accepted density for a maximum; 
#'                           as a multiple of kernel height
#' @return A list with min (minima), max (maxima), and
#'         detection.density (minimum detection density)
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples find_optima(rnorm(100), bw = 1)
#' @keywords utilities
find_optima <- function(f, detection = 0, bw = 1, detection.limit = 1) {

   # FIXME bw is now assumed to be 1. This may be far from
   # optimal. Should be determined automatically.

    # multiple of kernel height 
    kernel.height <- dnorm(0, sd = bw) / length(f) 
    deth <- detection * kernel.height 
    detl <- detection.limit * kernel.height 
    
    # Detect minima and maxima of the density (see Livina et al.) these correspond
    # to maxima and minima of the potential, respectively including end points of the
    # vector
    maxima <- find_maxima(f)
    minima <- find_minima(f)

    # remove maxima that are below detection limit
    maxima <- maxima[f[maxima] >= detl]
    minima <- remove_obsolete_minima(f, maxima, minima)
    minima <- unlist(minima)
    maxima <- unlist(maxima)

    # Remove minima and maxima that are too shallow
    delmini <- logical(length(minima))
    delmaxi <- logical(length(maxima))
    if (length(maxima) > 0) {
      for (j in 1:length(maxima)) {
        
        # Calculate distance of this maximum to all minima
        s <- minima - maxima[[j]]
        
        # Set distances to deleted minima to zero
        s[delmini] <- 0
        
        # identify the closest remaining minima
        i1 <- i2 <- NULL
        if (length(s) > 0) {
            
            minima.spos <- minima[s > 0]
            minima.sneg <- minima[s < 0]
            
            if (length(minima.spos) > 0) {
                i1 <- min(minima.spos)
            }
            if (length(minima.sneg) > 0) {
                i2 <- max(minima.sneg)
            }            
        }
        
        # if no positive differences available, set it to same value with i2
        if ((is.null(i1) && !is.null(i2))) {
            i1 <- i2
        } else if ((is.null(i2) && !is.null(i1))) {
            # if no negative differences available, set it to same value with i1
            i2 <- i1
        }
        
        if (!is.null(i1) && is.na(i1)) {
            i1 <- NULL
        }
        if (!is.null(i2) && is.na(i2)) {
            i2 <- NULL
        }
        
        # If a closest minimum exists, check differences and remove if difference is
        # under threshold
        if (!is.null(i1)) {
            
            # Smallest difference between this maximum and the closest minima
            diff <- min(c((f[maxima[[j]]] - f[i1]), (f[maxima[[j]]] - f[i2])))
            
            if (diff < deth) {
                
                # If difference is below threshold, delete this maximum
                delmaxi[[j]] <- TRUE
                
                # Delete the larger of the two neighboring minima
                if (f[[i1]] > f[[i2]]) {
                  delmini[minima == i1] <- TRUE
                } else {
                  delmini[minima == i2] <- TRUE
                }
            }
            
        } else {
            # if both i1 and i2 are NULL, do nothing
        }
      }
    } else {
      # Skip
    }
    # Delete the shallow minima and maxima
    if (length(minima) > 0 && sum(delmini) > 0) {
        minima <- minima[!delmini]
    }

    # Combine maxima that do not have minima in between
    if (length(maxima) > 1) {
      maxima2 <- c()
      for (i in 1:(length(maxima) - 1)) {
        nominima <- TRUE
	cnt <- 0
	while (nominima & (i + cnt) < length(maxima)) {
	cnt <- cnt + 1
	nominima <- sum(minima > maxima[[i]] & minima < maxima[[i + cnt]]) == 0
	# if (is.na(nominima)) {nominima <- TRUE}
      }
      maxs <- maxima[i:(i + cnt - 1)]
      maxima2 <- c(maxima2, round(mean(maxs[which(f[maxs] == max(f[maxs]))])))
    }
    if (!maxima[[length(maxima)]] %in% maxima2) {
    maxima2 <- c(maxima2, maxima[[length(maxima)]])
    }
    maxima <- maxima2
    }
   
    
    if (length(maxima) > 0 && sum(delmaxi) > 0) {
        maxima <- maxima[!delmaxi]
    }
    
    list(min = minima, max = maxima, detection.density = deth)
    
}

remove_obsolete_minima <- function (f, maxima, minima) {

  # remove minima that now became obsolete If there are multiple
  # minima between two consecutive maxima after removing the maxima
  # that did not pass the threshold, take the average of the minima;
  # return the list of indices such that between each pair of
  # consecutive maxima, there is exactly one minimum
  
    if (length(maxima) > 1) {
        minima <- sapply(2:length(maxima), function(i) {
            
            mins <- minima[minima >= maxima[[i - 1]] & minima <= maxima[[i]]]
            if (length(mins) > 0) {
                round(mean(mins[which(f[mins] == min(f[mins]))]))
            } else {
                NULL
            }
        })
        
    } else {
        minima <- NULL
    }

    # Remove minima that are outside the most extreme maxima
    minima <- minima[minima > min(maxima) & minima < max(maxima)]

    minima
}


 
find_minima <- function (f) {
  find_maxima(-f)
}

find_maxima <- function (f) {

    f2 <- c(Inf, -f, Inf)
    cnt <- 1
    ops <- c()
    opcnt <- 0
    while (cnt < length(f2)) {
      if (f2[[cnt + 1]] - f2[[cnt]] <= 0) {
	while (f2[[cnt + 1]] - f2[[cnt]] <= 0) {
	  cnt <- cnt + 1
        }
	ind1 <- cnt - 1
	while (f2[[cnt + 1]] - f2[[cnt]] == 0) {
	  cnt <- cnt + 1
        }
	if (f2[[cnt + 1]] - f2[[cnt]] > 0) {
	  ind2 <- cnt - 1
    	  opcnt <- opcnt + 1
	  ops[[opcnt]] <- round(mean(c(ind1, ind2)))
	} else if (f2[[cnt + 1]] - f2[[cnt]] < 0) {
	  ind2 <- NULL
	}
      }
      cnt <- cnt + 1
    }
    ops 
}




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

