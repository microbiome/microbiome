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
