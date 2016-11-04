#' @title Plot Potential
#' @description Visualization of the potential function.
#' @param res output from potential_slidingaverage function
#' @param cutoff parameter determining the upper limit of potential for
#'               visualizations
#' @param plot.contours Plot contour lines.
#' @param binwidth binwidth for contour plot
#' @param bins bins for contour plot. Overrides binwidth if given
#' @return A ggplot2 visualization of the potential landscape.
#' @export
#' @references Lahti et al. Tipping elements of the human intestinal ecosystem.
#'             \emph{Nature Communications} 5:4344, 2014. 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @details Applied on the output of the
#'          \code{\link{potential_slidingaverage}} function.
#' @examples
#'   X <- c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), 
#'   	    rnorm(1000, mean = 2))
#'   param <- seq(0,5,length=3000); 
#'   res <- potential_slidingaverage(X, param); 
#'   plot_potential(res$res, title = '', 
#'	       	    xlab.text = '', ylab.text = '', 
#'		    cutoff = 0.5, plot.contours = TRUE,
#'		    binwidth = 0.2)
#' @keywords utils
plot_potential <- function(res,
			   cutoff = 0.5,
	       	           plot.contours = TRUE, 
    			   binwidth = 0.2,
			   bins = NULL) {
    
    scale_fill_gradient <- NULL # Avoid build warnings
    
    cut.potential <- max(apply(res$pots, 1, min)) +
    		       cutoff * abs(max(apply(res$pots, 
        1, min)))  # Ensure all minima are visualized
    pots <- res$pots
    pots[pots > cut.potential] <- cut.potential
    
    # Static contour Interpolate potential grid
    intp <- tgp::interp.loess(as.vector(res$pars),
		              as.vector(res$xis), as.vector(pots), 
        gridlen = 2 * dim(pots))
    xy <- expand.grid(intp$x, intp$y)
    z <- as.vector(intp$z)
    z[is.na(z)] <- max(na.omit(z))
    
    potential <- NULL
    df <- data.frame(list(bg.var = xy[, 1], phylotype = xy[, 2], potential = z))
    bg.var <- NULL
    phylotype <- NULL
    
    p <- ggplot(df)
    p <- p + geom_tile(aes(x = bg.var, y = phylotype, z = potential,
      	                   fill = potential)) 
    p <- p + scale_fill_gradient(low = "black", high = "white")

    if (plot.contours) {
        if (!is.null(bins)) {
            warning("bins argument is overriding the binwidth argument!")
            p <- p + stat_contour(bins = bins,
	      	       aes(x = bg.var, y = phylotype, z = potential,
		       	     fill = potential))
        } else {
            p <- p + stat_contour(binwidth = binwidth,
	      	       aes(x = bg.var, y = phylotype, z = potential,
		       fill = potential))
        }
    }
    
    p
    
}

