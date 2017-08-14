#' @title Density Plot
#' @description Density visualization for data points overlaid on cross-plot.
#' @param x Data matrix to plot. The first two columns will be visualized as a cross-plot.
#' @param main title text
#' @param x.ticks Number of ticks on the X axis
#' @param rounding Rounding for X axis tick values
#' @param add.points Plot the data points as well
#' @param col Color of the data points. NAs are marked with darkgray.
#' @param adjust Kernel width adjustment
#' @param size point size
#' @param legend plot legend TRUE/FALSE
#' @return ggplot2 object
#' @examples p <- densityplot(cbind(rnorm(100), rnorm(100)))
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
densityplot <- function(x,
	       main = NULL,
	       x.ticks = 10,
	       rounding = 0, 
               add.points = TRUE,
	       col = "black",
	       adjust = 1,
	       size = 1,
	       legend = FALSE) {

    df <- x
    if (!is.data.frame(df)) {
      df <- as.data.frame(as.matrix(df))
    }

    # Avoid warnings
    x <- y <- ..density.. <- color <- NULL

    # If colors are NA:
    col <- as.character(col)
    col[unname(which(is.na(col)))] <- "darkgray"

    theme_set(theme_bw(20))
    xvar <- colnames(df)[[1]]
    yvar <- colnames(df)[[2]]
    df[["x"]] <- df[, 1]
    df[["y"]] <- df[, 2]
    df[["color"]] <- col
    df[["size"]] <- size

    # Remove NAs
    df <- df[!(is.na(df[["x"]]) | is.na(df[["y"]])), ]
    
    # Determine bandwidth for density estimation
    bw <- adjust * c(bandwidth.nrd(df[["x"]]), bandwidth.nrd(df[["y"]]))
    if (any(bw == 0)) {
      warning("Zero bandwidths (possibly due to small number of observations). Using minimal bandwidth.")
      bw[bw == 0] = bw[bw == 0] + min(bw[!bw == 0])
    }

    # Construct the figure
    p <- ggplot(df)
    p <- p + stat_density2d(aes(x, y, fill = ..density..), geom = "raster", h = bw, contour = FALSE)
    p <- p + scale_fill_gradient(low = "white", high = "black")


    if (add.points) {
      if (length(unique(df$color)) == 1 && length(unique(df$size)) == 1) {

        p <- p + geom_point(aes(x = x, y = y), col = unique(df$color), size = unique(df$size))
      } else if (length(unique(df$color)) == 1 && length(unique(df$size)) > 1) {
        p <- p + geom_point(aes(x = x, y = y, size = size), col = unique(df$color))
      } else if (length(unique(df$color)) > 1 && length(unique(df$size)) == 1) {
        p <- p + geom_point(aes(x = x, y = y, col = color), size = unique(df$size))
      } else {
        p <- p + geom_point(aes(x = x, y = y, col = color, size = size))
      }      
    }

    p <- p + xlab(xvar) + ylab(yvar)
    
    if (!legend) {
        p <- p + theme(legend.position = "none")
    }


    p <- p + scale_x_continuous(breaks = round(seq(floor(min(df[["x"]])), 
                   ceiling(max(df[["x"]])), 
                   length = x.ticks), rounding))

    if (!is.null(main)) {
        p <- p + ggtitle(main)
    }

    p
    
}
