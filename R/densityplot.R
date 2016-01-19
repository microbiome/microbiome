#' @title densityplot
#' @description Plots densities of data points in addition to cross-plot points.
#'
#' @param mat Data matrix to plot. The first two columns will be visualized 
#'              as a cross-plot.
#' @param main title text
#' @param x.ticks Number of ticks on the X axis
#' @param rounding Rounding for X axis tick values
#' @param add.points Plot the data points as well
#' @param col Color of the data points
#' @param adjust Kernel width adjustment
#' @param size point size
#' @param legend plot legend TRUE/FALSE
#'
#' @return ggplot2 object
#'
#' @examples p <- densityplot(cbind(rnorm(100), rnorm(100)))
#'
#' @export
#' @import ggplot2
#' @importFrom MASS bandwidth.nrd
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
densityplot <- function(mat, main = NULL, x.ticks = 10, rounding = 0, 
               add.points = TRUE, 
    col = "red", adjust = 1, size = 1, legend = FALSE) {

    mat <- as.matrix(mat)

    # Avoid warnings
    x <- y <- ..density.. <- color <- NULL
    
    theme_set(theme_bw(20))
    df <- as.data.frame(mat)
    xvar <- colnames(mat)[[1]]
    yvar <- colnames(mat)[[2]]
    df[["x"]] <- df[, 1]
    df[["y"]] <- df[, 2]
    df[["color"]] <- col
    df[["size"]] <- size
    
    # Remove NAs
    df <- df[!(is.na(df[["x"]]) | is.na(df[["y"]])), ]
    
    # Determine bandwidth for density estimation
    bw <- adjust * c(bandwidth.nrd(df[["x"]]), bandwidth.nrd(df[["y"]]))
    
    # Construct the figure
    p <- ggplot(df)
    p <- p + stat_density2d(aes(x, y, fill = ..density..), geom = "raster", h = bw, contour = FALSE)
    p <- p + scale_fill_gradient(low = "white", high = "black")
    
    if (add.points) {
        p <- p + geom_point(aes(x = x, y = y, col = color, size = size))
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
