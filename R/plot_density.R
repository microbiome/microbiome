#' @title Density plot
#' @description Plot abundance density across samples for a given taxon
#' @param x \code{\link{phyloseq-class}} object or an OTU matrix (samples x phylotypes)
#' @param variable OTU or metadata variable to visualize
#' @param log10 Logical. Show log10 abundances or not.
#' @param adjust see stat_density
#' @param kernel see stat_density
#' @param trim see stat_density
#' @param na.rm see stat_density
#' @param fill Fill color
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples # plot_density(x, variable = "Dialister")
#' @keywords utilities
plot_density <- function (x, variable = NULL, log10 = FALSE, adjust = 1, kernel = "gaussian", trim = FALSE, na.rm = FALSE, fill = "gray") {

  x <- pickdata(x, variable) 	     	   

  df <- data.frame(x = x)
  p <- ggplot(df, aes(x = x))
  p <- p + geom_density(adjust = adjust, kernel = kernel, trim = trim, na.rm = na.rm, fill = fill) 

  if (log10) {
    p <- p + scale_x_log10()  
    p <- p + ggtitle(paste(variable))
    p <- p + xlab("Abundance (Log10)")
  } else {
    p <- p + ggtitle(paste(variable))
    p <- p + xlab("Abundance (Abs.)")
  }

  p <- p + ylab("Frequency")

  p

}

