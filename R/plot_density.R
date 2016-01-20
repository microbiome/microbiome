#' @title Density plot
#' @description Plot abundance density across samples for a given taxon
#' @param x \code{\link{phyloseq-class}} object or an OTU matrix (samples x phylotypes)
#' @param otu.name OTU name to visualize
#' @param log10 Logical. Show log10 abundances or not.
#' @param adjust see stat_density
#' @param kernel see stat_density
#' @param trim see stat_density
#' @param na.rm see stat_density
#' @param fill Fill color
#' @return A \code{\link{ggplot}} plot object.
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @export
#' @examples # plot_density(x, otu.name = "Dialister")
#' @keywords utilities
plot_density <- function (x, otu.name = NULL, log10 = FALSE, adjust = 1, kernel = "gaussian", trim = FALSE, na.rm = FALSE, fill = "gray") {

  taxa_are_rows <- NULL

  if (is.vector(x)) { 

    df <- data.frame(x = x)
    p <- ggplot(df, aes(x = x))
    p <- p + geom_density(adjust = adjust, kernel = kernel, trim = trim, na.rm = na.rm, fill = fill)
    
    if (log10) {
      p <- p + scale_x_log10()
    }

  } else if (class(x) == "phyloseq") { 

    xx <- otu_table(x)@.Data

    if (!taxa_are_rows(x)) { xx <- t(xx)}
    
    xxx <- as.vector(xx[otu.name,])
    
    p <- plot_density(xxx, log10 = log10, adjust = adjust, kernel = kernel, trim = trim, na.rm = na.rm) 

  }

  if (log10) {
    p <- p + ggtitle(paste(otu.name))
    p <- p + xlab("Abundance (Log10)")
  } else {
    p <- p + ggtitle(paste(otu.name))
    p <- p + xlab("Abundance (Abs.)")
  }

  p <- p + ylab("Frequency")

  p

}

