#' @title Plot Density
#' @description Plot abundance density across samples for a given taxon.
#' @param x \code{\link{phyloseq-class}} object or an OTU matrix
#' (samples x phylotypes)
#' @param variable OTU or metadata variable to visualize
#' @param log10 Logical. Show log10 abundances or not.
#' @param adjust see stat_density
#' @param kernel see stat_density
#' @param trim see stat_density
#' @param na.rm see stat_density
#' @param tipping.point Optional. Indicate critical point for abundance
#' variations to be highlighted.
#' @param fill Fill color
#' @param xlim X axis limits
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples
#' # Load gut microbiota data on 1006 western adults
#' # (see help(atlas1006) for references and details)
#' data(dietswap)
#' # Use compositional abundances instead of absolute signal
#' pseq.rel <- transform(dietswap, 'compositional')
#' # Population density for Dialister spp.; with log10 on the abundance (X)
#' # axis
#' library(ggplot2)
#' p <- plot_density(pseq.rel, variable='Dialister') + scale_x_log10()
#' @keywords utilities
plot_density <- function(x, variable=NULL, log10=FALSE, adjust=1,
    kernel="gaussian", 
    trim=FALSE, na.rm=FALSE, fill="gray", tipping.point=NULL,
    xlim=NULL) {
    
    x <- pickdata(x, variable)
    
    if (log10 && min(x) == 0) {
        warning("The minimum abundance is 0. To avoid singularities with 
        log10, the abundances are shifted up by a small constant 
        (minimum non-zero value).")
        x=x + min(x[!is.na(x) & (x > 0)])
    }
    
    df <- data.frame(x=x)
    p <- ggplot(df, aes(x=x))
    
    p <- p + geom_density(adjust=adjust, kernel=kernel, trim=trim,
        na.rm=na.rm, 
        fill=fill)
    
    if (log10) {
        p <- p + scale_x_log10()
        p <- p + ggtitle(paste(variable))
        p <- p + labs(x = "Abundance (Log10)")
    } else {
        p <- p + ggtitle(paste(variable))
        p <- p + labs(x = "Abundance")
    }
    
    p <- p + labs(y = "Frequency")
    
    if (!is.null(tipping.point)) {
        p <- p + geom_vline(aes(xintercept=tipping.point), linetype=2,
                    size=1)
    }
    
    if (!is.null(xlim)) {
        p <- p + coord_cartesian(xlim)
    }
    
    p
    
}



pickdata <- function(x, otu.name) {
    
    if (is.vector(x)) {
        
        xxx <- x
        
    } else if (is.phyloseq(x)) {
        
        xx <- abundances(x)
        meta <- sample_data(x)
        
        # If OTU name not in otu data then try metadata
        if (otu.name %in% rownames(xx)) {
    
            xxx <- as.vector(xx[otu.name, ])
        
        } else if (otu.name %in% colnames(meta)) {
    
            xxx <- unlist(meta[, otu.name])
        
        } else {
    
            stop(paste("The provided variable name", otu.name,
                "cannot be found from OTUs or metadata"))
        
        }
    
    } else {

        stop("Provide proper variable for pickdata function.")
    
    }
    
    xxx
}
