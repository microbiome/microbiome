#' @title Bootstrapped Potential Analysis 
#' @description Analysis of multimodality based on bootstrapped potential
#' analysis of Livina et al. (2010) as described in Lahti et al. (2014).
#' @param x Input data vector
#' @param peak.threshold Mode detection threshold
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iter Bootstrap iterations
#' @param min.density minimum accepted density for a maximum; as a multiple of
#' kernel height
#' @return List with following elements:
#' \itemize{
#' \item{modes}{Number of modes for the input data vector
#' (the most frequent number of modes from bootstrap)}
#' \item{minima}{Average of potential minima across the bootstrap samples
#' (for the most frequent number of modes)}
#' \item{maxima}{Average of potential maxima across the bootstrap samples
#' (for the most frequent number of modes)}
#' \item{unimodality.support}{Fraction of bootstrap samples exhibiting
#' unimodality}
#' \item{bws}{Bandwidths}
#' }
#' @export
#' @examples
#' 
#' # Example data; see help(peerj32) for details
#' data(peerj32)
#' 
#' # Log10 abundance of Dialister
#' x <- abundances(transform(peerj32$phyloseq, "clr"))['Dialister',]
#'
#' # Bootstrapped potential analysis
#' # In practice, use more bootstrap iterations
#' # res <- potential_analysis(x, peak.threshold=0, bw.adjust=1,
#' #    bs.iter=9, min.density=1)
#'
#' @seealso plot_potential
#' @references
#' \itemize{
#' \item{}{Livina et al. (2010). Potential analysis reveals changing number
#' of climate states during the last 60 kyr.
#' \emph{Climate of the Past}, 6, 77-82.}
#' \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.}
#' }
potential_analysis <- function(x, peak.threshold=0, bw.adjust=1,
    bs.iter=100, min.density=1) {

    if (is.matrix(x) && nrow(x) == 1) {
        x <- as.vector(x)
    }
    
    # TODO :
    # Use densEstBayes() for density estimation and uncertainty analysis
    # https://arxiv.org/abs/2009.06182

    nmodes <- c()
    minpoints <- list()
    maxpoints <- list()
    bws <- c()
    
    for (r in seq_len(bs.iter)) {
        
        # Bootstrap
        rs <- sample(length(x), replace=TRUE)
        
        xbs <- na.omit(unname(x[rs]))
        
        a <- potential_univariate(xbs,
            grid.size=floor(0.2 * length(x)),
            peak.threshold=peak.threshold, 
            bw.adjust=bw.adjust, min.density=min.density)
        
        nmodes[[r]] <- length(a$max.points)
        minpoints[[r]] <- a$min.points
        maxpoints[[r]] <- a$max.points
        bws[[r]] <- a$bw
        
    }

    nmodes <- unlist(nmodes)

    # Most frequently observed number of modes
    top.modes  <- as.numeric(names(which.max(table(nmodes))))

    min.points <- colMeans(do.call("rbind", minpoints[nmodes == top.modes]))
    max.points <- colMeans(do.call("rbind", maxpoints[nmodes == top.modes]))

    unimodality.support <- mean(nmodes <= 1)

    # Return the most frequent number of modes and the corresponding
    # tipping points from the bootstrap analysis
    list(modes=top.modes,
        minima=min.points, maxima=max.points,
        unimodality.support=unimodality.support, 
        bws=bws)
    
}




#' @title Potential Analysis for Univariate Data
#' @description One-dimensional potential estimation for univariate timeseries.
#' @param x Univariate data (vector) for which the potentials shall be estimated
#' @param std Standard deviation of the noise (defaults to 1; this will set
#'    scaled potentials)
#' @param bw kernel bandwidth estimation method 
#' @param weights optional weights in ksdensity
#' (used by potential_slidingaverages).
#' @param grid.size Grid size for potential estimation.
#' of density kernel height dnorm(0, sd=bandwidth)/N
#' @param bw.adjust The real bandwidth will be bw.adjust*bw; defaults to 1
#' @param density.smoothing Add a small constant density across the
#' whole observation range to regularize density estimation (and to
#' avoid zero probabilities within the observation range). This
#' parameter adds uniform density across the observation range, scaled
#' by density.smoothing.
#' @param min.density minimum accepted density for a maximum; as a
#' multiple of kernel height
#' @inheritParams potential_analysis
#' @return \code{potential_univariate} returns a list with the
#' following elements:
#' \itemize{
#' \item{xi }{the grid of points on which the potential is estimated}
#' \item{pot }{The estimated potential: -log(f)*std^2/2,
#' where f is the density.}
#' \item{density }{Density estimate corresponding to the potential.}
#' \item{min.inds }{indices of the grid points at which the density has
#' minimum values; (-potentials; neglecting local optima)}
#' \item{max.inds }{indices the grid points at which the density has
#' maximum values; (-potentials; neglecting local optima)}
#' \item{bw }{bandwidth of kernel used}
#' \item{min.points }{grid point values at which the density has
#' minimum values; (-potentials; neglecting local optima)}
#' \item{max.points }{grid point values at which the density has
#' maximum values; (-potentials; neglecting local optima)}
#' }
#' @references
#' \itemize{
#' \item{}{Livina et al. (2010).
#'  Potential analysis reveals changing number of climate states during
#'  the last 60 kyr. \emph{Climate of the Past}, 6, 77-82.}
#' \item{}{Lahti et al. (2014).
#'  Tipping elements of the human intestinal ecosystem.
#'  \emph{Nature Communications} 5:4344.}
#' }
#' @author Based on Matlab code from Egbert van Nes modified by Leo Lahti.
#' Extended from the initial version in the \pkg{earlywarnings} R package.
#' @examples 
#' # res <- potential_univariate(x)
#' @keywords early-warning
potential_univariate <- function(x, std=1, bw="nrd", weights=c(),
    grid.size=NULL, 
    peak.threshold=1, bw.adjust=1, density.smoothing=0, min.density=1) {
    
    if (is.null(grid.size)) {
        grid.size <- floor(0.2 * length(x))
    }
    
    # Density estimation
    tmp <- try(de <- density(x, bw=bw, adjust=bw.adjust,
        kernel="gaussian", 
        weights=weights,
        window=kernel,
        n=grid.size,
        from=min(x), to=max(x), 
        cut=3, na.rm=FALSE))
    if (length(tmp) == 1 && is(tmp) == "try-error") {
        # Just use default parameters if failing otherwise
        warning("Density estimation with custom parameters failed. 
            Using the defaults.")
        de <- density(x)
    }
    
    # Smooth the estimated density (f <- de$y) by adding a small
    # probability across
    # the whole observation range (to avoid zero probabilities for points
    # in the observation range)
    f <- de$y + density.smoothing * 1/diff(range(de$x))  # *max(de$y)
    
    # Normalize the density such that it integrates to unity
    f <- f/sum(diff(de$x[seq_len(2)]) * f)
    
    # Final grid points and bandwidth
    grid.points <- de$x
    bw <- de$bw
    
    # Compute potential
    U <- -log(f) * std^2/2
    
    # Backtransform to density distribution f <- exp(-2*U/std^2)
    
    # Ignore very local optima
    
    # Note mins and maxs for density given # here (not for potential, which h
    # has the opposite signs)
    ops <- find_optima(f, peak.threshold=peak.threshold, bw=bw,
        min.density=min.density)
    min.points <- grid.points[ops$min]
    max.points <- grid.points[ops$max]
    peak.threshold2 <- ops$peak.threshold2
    
    list(grid.points=grid.points, pot=U, density=f, min.inds=ops$min,
        max.inds=ops$max, 
        bw=bw, min.points=min.points, max.points=max.points,
    peak.threshold2=peak.threshold2)
    
}





#' @title Find Optima
#' @description Detect optima, excluding local optima below peak.threshold. 
#' @param f density
#' @param bw bandwidth
#' @param min.density Minimun accepted density for a maximum; 
#'                        as a multiple of kernel height
#' @inheritParams potential_analysis
#' @return A list with min (minima), max (maxima), and
#'        peak.threshold (minimum detection density)
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#'    # Not exported
#'    # o <- find_optima(rnorm(100), bw=1)
#' @keywords utilities
find_optima <- function(f, peak.threshold=0, bw=1, min.density=1) {
    
    # FIXME bw is now assumed to be 1. This may be far from optimal. Should be
    # determined automatically.
    
    # multiple of kernel height
    kernel.height <- dnorm(0, sd=bw)/length(f)
    peak.threshold2 <- peak.threshold * kernel.height
    detl <- min.density * kernel.height
    
    # Detect minima and maxima of the density (see Livina et al.)
    # these correspond to
    # maxima and minima of the potential, respectively including end
    # points of the
    # vector
    maxima <- find_maxima(f)
    minima <- find_minima(f)
    
    # remove maxima that are below min.density
    maxima <- maxima[f[maxima] >= detl]

    minima <- remove_obsolete_minima(f, maxima, minima)
    minima <- unlist(minima)
    maxima <- unlist(maxima)
    
    # Remove minima and maxima that are too shallow
    delmini <- logical(length(minima))
    delmaxi <- logical(length(maxima))
    if (length(maxima) > 0) {
        for (j in seq_len(length(maxima))) {
            
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
            
            # if no positive differences available, set it to
        # same value with i2
            if ((is.null(i1) && !is.null(i2))) {
                i1 <- i2
            } else if ((is.null(i2) && !is.null(i1))) {
                # if no negative differences available,
        # set it to same value with i1
                i2 <- i1
            }
            
            if (!is.null(i1) && is.na(i1)) {
                i1 <- NULL
            }
            if (!is.null(i2) && is.na(i2)) {
                i2 <- NULL
            }
            
            # If a closest minimum exists, check differences and
        # remove if difference is
            # under threshold
            if (!is.null(i1)) {
                
                # Smallest difference between this maximum and the
                # closest minima
                diff <- min(c((f[maxima[[j]]] - f[i1]),
            (f[maxima[[j]]] - f[i2])))
                
                if (diff < peak.threshold2) {
                
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
        for (i in seq_len((length(maxima) - 1))) {
            nominima <- TRUE
            cnt <- 0
            while (nominima & (i + cnt) < length(maxima)) {
                cnt <- cnt + 1
                nominima <- sum(minima > maxima[[i]] &
            minima < maxima[[i + cnt]]) == 0
                
                # if (is.na(nominima)) {nominima <- TRUE}
            }
            maxs <- maxima[i:(i + cnt - 1)]
            maxima2 <- c(maxima2,
            round(mean(maxs[which(f[maxs] == max(f[maxs]))])))
        }
        if (!maxima[[length(maxima)]] %in% maxima2) {
            maxima2 <- c(maxima2, maxima[[length(maxima)]])
        }
        maxima <- maxima2
    }
    
    
    if (length(maxima) > 0 && sum(delmaxi) > 0) {
        maxima <- maxima[!delmaxi]
    }
    
    list(min=minima, max=maxima, peak.threshold2=peak.threshold2)
    
}

remove_obsolete_minima <- function(f, maxima, minima) {
    
    # remove minima that now became obsolete If there are multiple
    # minima between two consecutive maxima after removing the maxima
    # that did not pass the threshold, take the average of the minima;
    # return the list of indices such that between each pair of
    # consecutive maxima, there is exactly one minimum
    
    if (length(maxima) > 1) {
        minima <- vapply(2:length(maxima), function(i) {
            
            mins <- minima[minima >= maxima[[i - 1]] & minima <= maxima[[i]]]
            if (length(mins) > 0) {
                round(mean(mins[which(f[mins] == min(f[mins]))]))
            } else {
                NULL
            }
        }, 1)
        
    } else {
        minima <- NULL
    }
    
    # Remove minima that are outside the most extreme maxima
    minima <- minima[minima > min(maxima) & minima < max(maxima)]
    
    minima
}



find_minima <- function(f) {
    find_maxima(-f)
}

find_maxima <- function(f) {
    
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
    unlist(ops)
}




