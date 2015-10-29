#' @title plot_density
#' @description Plot taxon density for samples
#'
#' @param x \code{\link{phyloseq-class}} object or an OTU matrix (samples x phylotypes)
#' @param otu.name OTU name to visualize
#' @param log10 Logical. Show log10 abundances or not.
#' @param adjust see stat_density
#' @param kernel see stat_density
#' @param trim see stat_density
#' @param na.rm see stat_density
#' @param fill Fill color
#'
#' @return A \code{\link{ggplot}} plot object.
#' 
#' @import ggplot2
#' @export
#' @examples # 
#' @keywords utilities
plot_density <- function (x, otu.name = NULL, log10 = FALSE, adjust = 1, kernel = "gaussian", trim = FALSE, na.rm = FALSE, fill = "gray") {

  if (is.vector(x)) { 

    df <- data.frame(x = x)
    p <- ggplot(df, aes(x = x))
    p <- p + geom_density(adjust = adjust, kernel = kernel, trim = trim, na.rm = na.rm, fill = fill)
    
    if (log10) {
      #x <- log10(x)
      p <- p + scale_x_log10()
    }

  } else if (class(x) == "phyloseq") { 

    x <- otu_table(x)@.Data[otu.name,]
    p <- plot_density(x, log10 = log10, adjust = adjust, kernel = kernel, trim = trim, na.rm = na.rm) 

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

