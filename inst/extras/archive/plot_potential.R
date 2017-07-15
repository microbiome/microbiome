#' @title Plot Potential
#' @description Visualization of the potential function.
#' @param res output from potential_slidingaverage function
#' @param cutoff parameter determining the upper limit of potential for
#' visualizations
#' @param plot.contours Plot contour lines.
#' @param binwidth binwidth for contour plot
#' @param bins bins for contour plot. Overrides binwidth if given
#' @return A ggplot2 visualization of the potential landscape.
#' @export
#' @references Lahti et al. Tipping elements of the human intestinal ecosystem.
#' \emph{Nature Communications} 5:4344, 2014. 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @details Applied on the output of the
#' \code{\link{potential_slidingaverage}} function.
#' @examples
#' # Load gut microbiota data on 1006 western adults
#' # From http://doi.org/10.5061/dryad.pk75d
#' # (see help(atlas1006) for references and details)
#' data(atlas1006)
#' pseq <- subset_samples(atlas1006, DNA_extraction_method == 'r')
#'
#' # Pick diversity and age
#' diversity <- exp(global(pseq, index='shannon'))
#' age <- meta(pseq)$age
#'   
#' # Run potential analysis
#' # library(earlywarnings)
#' # library(ggplot2)
#' # suppressWarnings(res <- movpotential_ews(diversity, age)$res)
#' # suppressWarnings(plot_potential(res, cutoff=0.5))
#' @seealso potential_analysis
#' @keywords utils
plot_potential <- function(res, cutoff=0.5, plot.contours=TRUE,
    binwidth=0.2, bins=NULL) {
    
    
    # TODO Experimental implementation including bootstrap as in Lahti
    # et al (2014) Nat Comm potential_analysis(diversity, age) to
    # replace movpotential_ews, see potential_analysis.R does not work
    # with visualization since not meant for that purpose TODO enable
    # plotting from potential_analysis omitting bootstrap to get rid
    # of earlywarnings dependency
    
    scale_fill_gradient <- NULL  # Avoid build warnings
    
    potentials <- res$pots
    covariates <- res$pars
    grid <- res$xis
    
    cut.potential <- max(apply(potentials, 1, min)) +
        cutoff * abs(max(apply(potentials, 
        1, min)))  # Ensure all minima are visualized
    
    potentials[potentials > cut.potential] <- cut.potential
    
    # Static contour Interpolate potential grid
    intp <- loess_interpolation(as.vector(covariates),
        as.vector(grid), as.vector(potentials), 
        gridn=2 * dim(potentials))
    
    xy <- expand.grid(intp$x, intp$y)
    z <- as.vector(intp$z)
    z[is.na(z)] <- max(na.omit(z))
    
    potential <- NULL
    df <- data.frame(list(bg.var=xy[, 1], phylotype=xy[, 2], potential=z))
    bg.var <- NULL
    phylotype <- NULL
    
    p <- ggplot(df)
    p <- p + geom_tile(aes(x=bg.var, y=phylotype,
        z=potential, fill=potential))
    p <- p + scale_fill_gradient(low="black", high="white")
    
    if (plot.contours) {
        if (!is.null(bins)) {
            warning("bins argument is overriding the binwidth argument!")
            p <- p + stat_contour(bins=bins,
            aes(x=bg.var, y=phylotype, z=potential, 
                fill=potential))
        } else {
            p <- p + stat_contour(binwidth=binwidth,
            aes(x=bg.var, y=phylotype, 
                z=potential, fill=potential))
        }
    }
    
    p
    
}




loess_interpolation <- function(x, y, z, gridn=c(40, 40), span=0.1, ...) {
    
    # Inspired by tgp::interp.loess but rewritten
    
    if (length(gridn) == 1) 
        gridn <- rep(gridn, 2)
    if (length(gridn) != 2) 
        stop("gridn is not 2")
    if (length(x) != length(y) && length(y) != length(z)) 
        stop("x, y, z lengths are not equal")
    if (length(x) < 30 && span < 0.5) {
        warning("with less than 30 points, 
        use span >> 0.5 or use akima", immediate.=TRUE)
        message(paste("trying span ", span, "for", length(x), "points\n"))
    }
    x0 <- seq(min(x), max(x), length=gridn[[1]])
    y0 <- seq(min(y), max(y), length=gridn[[2]])
    
    coords <- suppressWarnings(loess(z ~ x + y, data.frame(x=x, y=y),
        span=span, 
        ...))
    g <- expand.grid(x=x0, y=y0)
    g.pred <- predict(coords, g)
    
    list(x=x0, y=y0, z=matrix(g.pred, nrow=gridn))
    
}


#' @title Moving Average Potential
#' @description This function reconstructs a potential derived from data
#'            along a gradient of a given parameter.
#' @param X a vector of the X observations of the state variable of interest
#' @param param parameter values corresponding to the observations in X 
#' @param bw Bandwidth for smoothing kernels. Automatically determined
#'            by default.
#' @param bw.adjust Bandwidth adjustment constant
#' @param std Standard deviation.
#' @param grid.size number of evaluation points; number of steps between
#'         min and max potential; also used as kernel window size
#' @param plot.cutoff cuttoff for potential minima and maxima in visualization
#' @param plot.contours Plot contours on the landscape visualization
#' @param binwidth binwidth for contour plot
#' @param bins bins for contour plot. Overrides binwidth if given
#' @inheritParams potential_analysis
#' @return A list with the following elements:
#' \itemize{
#'    \item{pars}{values of the covariate parameter as matrix}
#'    \item{xis}{values of the x as matrix}
#'    \item{pots}{smoothed potentials}
#'    \item{mins}{minima in the densities
#'        (-potentials; neglecting local optima)}
#'    \item{maxs}{maxima in densities (-potentials; neglecting local optima)}
#'    \item{plot}{an object that displays the potential estimated in 2D}
#' }
#' @references
#'  \itemize{
#'   \item{}{Hirota, M., Holmgren, M., van Nes, E.H. & Scheffer, M. (2011).
#'            Global resilience of tropical forest and savanna to critical
#'            transitions. \emph{Science}, 334, 232-235.}
#'   \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#'        ecosystem. \emph{Nature Communications} 5:4344.}
#'  }
#' @author Leo Lahti, adapted from original Matlab code by Egbert van Nes.
#' @seealso \code{potential_univariate}
#' @examples
#'   \dontrun{
#'    # Not exported
#'    X <- c(rnorm(1000, mean=0),
#'            rnorm(1000, mean=-2),
#'            rnorm(1000, mean=2));
#'        param=seq(0,5,length=3000); 
#'        res <- potential_slidingaverage(X, param)
#'    }
#' @keywords utils
potential_slidingaverage <- function(X, param=NULL, bw="nrd",
    bw.adjust=1, 
    peak.threshold=0.1, std=1, grid.size=50, plot.cutoff=0.5,
    plot.contours=TRUE, 
    binwidth=0.2, bins=NULL) {
    
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
    xi <- seq(min(X), max(X), length=grid.size)
    
    # Initialize
    xis <- pars <- pots <- matrix(0, nrow=grid.size, ncol=length(xi))
    maxs <- mins <- matrix(0, nrow=grid.size, ncol=length(xi))
    
    for (i in 1:grid.size) {
        
        # Increase the parameter at each step
        par <- minparam + (i - 0.5) * step
        
        # Check which elements in evaluation range (param) are within 2*sd
        # of par
        weights <- exp(-0.5 * (abs(par - param)/sdwindow)^2)
        
        # LL: Normalization was added in the R implementation 16.5.2012
        weights <- weights/sum(weights)
        
        # Calculate the potential
        tmp <- potential_univariate(x=X, std=std, bw=bw,
        bw.adjust=bw.adjust, 
            weights=weights, grid.size=grid.size)
        
        # Store variables
        pots[i, ] <- tmp$pot
        xis[i, ] <- tmp$grid.points
        pars[i, ] <- par + rep(0, length(tmp$grid.points))
        mins[i, tmp$min.inds] <- 1
        maxs[i, tmp$max.inds] <- 1
        
    }
    
    res <- list(pars=pars, xis=xis, pots=pots, mins=mins,
        maxs=maxs, std=std)
    
    p <- plot_potential(res, cutoff=plot.cutoff,
        plot.contours=plot.contours, 
        binwidth=binwidth, bins=bins)
    
    p <- p + xlab("parameter/time") + ylab("state variable")
    
    list(res=res, plot=p)
    
}

