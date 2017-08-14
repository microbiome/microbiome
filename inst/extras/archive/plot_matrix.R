#' @title Matrix Heatmap
#' @description Fast investigation of matrix objects;
#' standard visualization choices made automatically.
#' @param mat matrix
#' @param type String. Specifies visualization type. Options: 'oneway' 
#' (color scale ranges from white to dark red; 
#' the color can be changed if needed); 
#' 'twoway' (color scale ranges from dark blue 
#' through white to dark red; colors can be changed 
#' if needed)
#' @param midpoint middle point for the color plot: smaller values are 
#' shown with blue, larger are shown with red 
#' in type='twoway'
#' @param palette Optional. Color palette.
#' @param colors Optional. Colors.
#' @param col.breaks breakpoints for the color palette
#' @param interval interval for palette color switches
#' @param plot_axes String. Indicates whether to plot 
#' x-axis ('x'), y-axis ('y'), or both ('both').
#' @param row.tick interval for plotting row axis texts
#' @param col.tick interval for plotting column axis texts
#' @param cex.xlab use this to specify distinct font size for the x axis
#' @param cex.ylab use this to specify distinct font size for the y axis
#' @param xlab optional x axis labels
#' @param ylab optional y axis labels
#' @param limit.trunc color scale rounding
#' @param cap Color scale end point
#' @param mar image margins
#' @param ... optional parameters to be passed to function 'image',
#' see help(image) for further details
#' @return A list with the color palette (colors), 
#' color breakpoints (breaks), and palette function (palette.function)
#' @export
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples 
#' mat <- rbind(c(1,2,3,4,5), c(1, 3, 1), c(4,2,2))
#' res <- plot_matrix(mat, 'twoway', midpoint=3) 
#' @keywords utilities
plot_matrix <- function(mat, type="twoway", midpoint=0, palette=NULL,
    colors=NULL, col.breaks=NULL, interval=0.1, plot_axes="both",
    row.tick=1, col.tick=1, cex.xlab=0.9, cex.ylab=0.9, xlab=NULL,
    ylab=NULL, limit.trunc=0, cap=NULL, mar=c(5, 4, 4, 2), ...) {
    
    # Center the data and color breakpoints around the specified midpoint
    mat <- mat - midpoint
    if (is.null(cap)) {
        cap <- max(abs(mat))
    }
    
    mat[mat > cap] <- cap
    mat[mat < -cap] <- (-cap)
        
    if (length(col.breaks) == 0) {
        
        m <- max(round(max(abs(mat)), limit.trunc) - interval, 0)
        
        mm <- m + interval/2
        
        vals <- seq(interval/2, mm, interval)
        
        # Set col.breaks evenly around zero
        col.breaks <- c(-(m + 1e+06), c(-rev(vals), vals), m + 1e+06)
        
    }
    
    my.palette <- palette
    if (is.null(palette)) {
        my.palette <- colorRampPalette(c("blue", "white", "red"), space="rgb")
    } else if (!class(palette) == "function" && palette == "blue-black-red") {
        my.palette <- colorRampPalette(c("blue", "black", "red"), space="rgb")
    } else if (!class(palette) == "function" && palette == "blue-white-red") {
        my.palette <- colorRampPalette(c("blue", "white", "red"), space="rgb")
    } else if (!class(palette) == "function" &&
        palette == "blue-white-yellow") {
        my.palette <- colorRampPalette(c("blue", "white", "yellow"),
    space="rgb")
    } else if (!class(palette) == "function" &&
        palette == "blue-black-yellow") {
        my.palette <- colorRampPalette(c("blue", "black", "yellow"),
    space="rgb")
    } else if (!class(palette) == "function" && palette == "bw") {
        gray.palette <- function(int) {
            gray(seq(0, 1, length=int))
        }
        my.palette <- gray.palette
    }
    
    palette <- my.palette
    
    # if mycolors is provided it overrides palette
    if (is.null(colors)) {
        colors <- palette(length(col.breaks) - 1)
    }
    
    # transpose and revert row order to plot matrix in the same way it
    # appears in its numeric form
    
    par(mar=mar)
    
    matm <- matrix(mat[rev(seq(nrow(mat))), ], ncol=ncol(mat))
    dimnames(matm) <- dimnames(mat)
    mat <- matm
    
    nsamples <- ncol(mat)
    nfeats <- nrow(mat)
    if (nfeats == 1) {
        mat <- mat
    }
    image(t(mat), col=colors, xaxt="n", yaxt="n",
        zlim=range(col.breaks), 
        breaks=col.breaks, ...)
    
    if (plot_axes == "both" || plot_axes == TRUE) {
        
        if (is.null(xlab)) {
            v <- seq(1, nsamples, col.tick)  # take every nth index
            axis(1, at=seq(0, 1, length=nsamples)[v],
        labels=colnames(mat)[v], 
            cex.axis=cex.xlab, las=2, ...)
        } else if (!xlab == "") {
            axis(1, at=seq(0, 1, length=length(xlab)),
            labels=xlab, cex.axis=cex.xlab, 
                las=2, ...)
        }
        
        if (is.null(ylab)) {
            v <- seq(1, nfeats, row.tick)  # take every nth index
            axis(2, at=seq(0, 1, length=nfeats)[v],
            labels=rev(rownames(mat))[v], 
                cex.axis=cex.ylab, las=2, ...)
        } else if (!ylab == "") {
            axis(2, at=seq(0, 1, length=nfeats),
            labels=rev(ylab), cex.axis=cex.ylab, 
                las=2, ...)
        }
        
    } else if (plot_axes == "x") {
        
        if (is.null(xlab)) {
            v <- seq(1, nsamples, col.tick)  # take every nth index
            axis(1, at=seq(0, 1, length=nsamples)[v],
            labels=colnames(mat)[v], 
                cex.axis=cex.ylab, las=2)
        } else {
            axis(1, at=seq(0, 1, length=nsamples),
            labels=ylab, cex.axis=cex.ylab, 
                las=2)
        }
        
    } else if (plot_axes == "y") {
        
        if (is.null(ylab)) {
            v <- seq(1, nfeats, row.tick)  # take every nth index
            axis(2, at=seq(0, 1, length=nfeats)[v],
            labels=rev(rownames(mat))[v], 
                cex.axis=cex.xlab, las=2)
        } else {
            axis(2, at=seq(0, 1, length=nfeats),
            labels=ylab, cex.axis=cex.xlab, 
                las=2)
        }
    }
    
    # Return default margins
    par(mar=c(5, 4, 4, 2) + 0.1)
    return(list(colors=colors, breaks=col.breaks +
        midpoint, palette.function=palette))
    
}

