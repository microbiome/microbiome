#' @title Find Optima
#' @description Detect optima, excluding local optima below
#'              detection.threshold. 
#' @param f density
#' @param detection.threshold detection threshold for peaks
#' @param bw bandwidth
#' @param detection.limit Minimun accepted density for a maximum; 
#'                           as a multiple of kernel height
#' @return A list with min (minima), max (maxima), and
#'         detection.density (minimum detection density)
#' @export
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples find_optima(rnorm(100), bw = 1)
#' @keywords utilities
find_optima <- function(f, detection.threshold = 0, bw = 1, detection.limit = 1) {

   # FIXME bw is now assumed to be 1. This may be far from
   # optimal. Should be determined automatically.

    # multiple of kernel height 
    kernel.height <- dnorm(0, sd = bw) / length(f) 
    deth <- detection.threshold * kernel.height 
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


