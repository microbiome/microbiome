# "A model is a lie that helps you see the truth."
#                                 - Howard Skipper

# Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' DensityPlot
#'
#' Description: Plots densities of data points in addition to cross-plot points.
#'
#' Arguments:
#'   @param mat Data matrix to plot. The first two columns will be visualized as a cross-plot.
#'   @param main title text
#'   @param x.ticks Number of ticks on the X axis
#'   @param rounding Rounding for X axis tick values
#'   @param add.points Plot the data points as well
#'   @param col Color of the data points
#'   @param adjust Kernel width adjustment
#'   @param size point size
#'
#' Returns:
#'   @return ggplot2 object
#'
#' @examples p <- densityplot(cbind(rnorm(100), rnorm(100)))
#'
#' @export
#' @import ggplot2
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

densityplot <- function (mat, main = NULL, x.ticks = 10, rounding = 0, add.points = TRUE, col = "red", adjust = 1, size = 1) {

    # mat: samples x features data matrix	     

    # Avoid warnings
    x <- y <- ..density.. <- color <- NULL

    theme_set(theme_bw(20))
    df <- as.data.frame(mat)
    xvar <- colnames(mat)[[1]]
    yvar <- colnames(mat)[[2]]
    df[["x"]] <- df[, 1]		
    df[["y"]] <- df[, 2]		
    df[["color"]] <- col

    # Remove NAs
    df <- df[!(is.na(df[["x"]]) | is.na(df[["y"]])), ]

    # Determine bandwidth for density estimation
    InstallMarginal("MASS")
    bw <- adjust*c(bandwidth.nrd(df[["x"]]), bandwidth.nrd(df[["y"]]))

    # Construct the figure
    p <- ggplot(df) 
    p <- p + stat_density2d(aes(x, y, fill=..density..), geom="raster", stat_params = list(h = bw, contour = F), geom_params = list()) 
    p <- p + scale_fill_gradient(low="white", high="black") 

    # MASS::kde2d(x, y)
    #des <- MASS::kde2d(x, y)
    #x <- des$x
    #y <- des$y
    #z <- des$z
    #df <- data.frame(cbind(expand.grid(des$x, des$y), as.vector(des$z)))
    #colnames(df) <- c("x", "y", "z")
    #p <- ggplot(df) 
    #p <- p + stat_density2d(aes(x, y, fill = z), geom = "raster", stat_params = list(h = bw, contour = F), geom_params = list()) 
    #p <- p + scale_fill_gradient(low = "white", high = "black") 

    if (add.points) {
      p <- p + ggplot2::geom_point(aes(x = x, y = y, col = color), size = size) 
    }

    p <- p + ggplot2::xlab(xvar) + ggplot2::ylab(yvar) 

    p <- p + ggplot2::theme(legend.position="none")
    p <- p + ggplot2::scale_x_continuous(breaks = round(seq(floor(min(df[["x"]])), ceiling(max(df[["x"]])), length = x.ticks), rounding))

    if (!is.null(main)) {
      p <- p + ggplot2::ggtitle(main)
    }

    p

}



#' correlation.heatmap
#'
#' Description: Visualizes n x m correlation table as heatmap. See examples for details.
#'
#' Arguments:
#'   @param df Data frame. Each row corresponds to a pair of correlated variables. The columns give variable names, correlations and significance estimates.
#'   @param Xvar X axis variable column name. For instance "X".
#'   @param Yvar Y axis variable column name. For instance "Y".
#'   @param fill Column to be used for heatmap coloring. For instance "correlation".
#'   @param star Column to be used for cell highlighting. For instance "p.adj".
#'   @param p.adj.threshold Significance threshold for the stars.
#'   @param correlation.threshold Include only elements that have absolute correlation higher than this value
#'   @param step color interval
#'   @param colours heatmap colours
#'   @param limits colour scale limits
#'   @param legend.text legend text
#'   @param order.rows Order rows to enhance visualization interpretability
#'   @param order.cols Order columns to enhance visualization interpretability
#'   @param text.size Adjust text size
#'   @param filter.significant Keep only the elements with at least one significant entry
#'
#' Returns:
#'   @return ggplot2 object
#'
#' @import ggplot2 
#' @importFrom ggplot2 theme_set
#'
#' @examples data(peerj32); cc <- cross.correlate(peerj32$lipids[, 1:10], peerj32$microbes[, 1:10]); p <- correlation.heatmap(cc, "X1", "X2", "Correlation")
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

correlation.heatmap <- function (df, Xvar, Yvar, fill, star = "p.adj",
p.adj.threshold = 1, correlation.threshold = 0, step = 0.2, colours =
c("darkblue", "blue", "white", "red", "darkred"), limits = NULL,
legend.text = "", order.rows = TRUE, order.cols = TRUE, text.size = 10, filter.significant = TRUE) {

  # df <- cc; Xvar <- "X1"; Yvar <- "X2"; fill = "cor"; star = "p.adj"; p.adj.threshold = 1; correlation.threshold = 0; order.rows = TRUE; order.cols = TRUE; text.size = 12; filter.significant = TRUE; step = 0.2; colours = c("darkblue", "blue", "white", "red", "darkred"); limits = c(-1, 1); legend.text = fill

  if (is.null(limits)) { 
    maxval <- max(abs(df[[fill]])) 
    if (maxval <= 1) {
      limits <- c(-1, 1)
    } else {
      limits <- c(-maxval, maxval)
    }
  }


  if (nrow(df) == 0) {warning("Input data frame is empty."); return(NULL)}

  if (filter.significant) {
    keep.X <- as.character(unique(df[((df[[star]] < p.adj.threshold) & (abs(df[[fill]]) > correlation.threshold)), Xvar]))
    keep.Y <- as.character(unique(df[((df[[star]] < p.adj.threshold) & (abs(df[[fill]]) > correlation.threshold)), Yvar]))
    df <- df[((df[[Xvar]] %in% keep.X) & (df[[Yvar]] %in% keep.Y)),]
  }		    

  theme_set(theme_bw(text.size))

  if (any(c("XXXX", "YYYY", "ffff") %in% names(df))) { stop("XXXX, YYYY, ffff are not allowed in df") }

  # ------------------------------------------------------

  df[[Xvar]] <- factor(df[[Xvar]])
  df[[Yvar]] <- factor(df[[Yvar]])

  if (order.rows || order.cols) {

    rnams <- unique(as.character(df[[Xvar]]))
    cnams <- unique(as.character(df[[Yvar]]))
 
    mat <- matrix(0, nrow = length(rnams), ncol = length(cnams))
    rownames(mat) <- rnams
    colnames(mat) <- cnams
    for (i in 1:nrow(df)) {
      mat[as.character(df[i, Xvar]), as.character(df[i, Yvar])] <- df[i, fill]
    }

    hm <- heatmap(mat)
    dev.off()
    rind <- hm$rowInd
    cind <- hm$colInd

    if (order.cols) {
      message("Ordering columns")
      df[[Xvar]] <- factor(df[[Xvar]], levels = rownames(mat)[rind])
    } 

    if (order.rows) {
      message("Ordering rows")
      df[[Yvar]] <- factor(df[[Yvar]], levels = colnames(mat)[cind])
    }

  }

  # ---------------------------------------------

  XXXX <- YYYY <- ffff <- NULL

  df[["XXXX"]] <- df[[Xvar]]
  df[["YYYY"]] <- df[[Yvar]]
  df[["ffff"]] <- df[[fill]]

  p <- ggplot(df, aes(x = XXXX, y = YYYY, fill = ffff))

  p <- p + geom_tile()

  p <- p + scale_fill_gradientn(legend.text, 
       	   		breaks = seq(from = min(limits), to = max(limits), by = step), 
  			colours = colours, 
  			limits = limits)

  p <- p + xlab("") + ylab("")

  p <- p + theme(axis.text.x = element_text(angle = 90))

  # Mark significant cells with stars
  inds <- which((df[[star]] < p.adj.threshold) & (abs(df[[fill]]) > correlation.threshold))
  if (!is.null(star) & length(inds) > 0) {
    df.sub <- df[inds,] 
    p <- p + geom_text(data = df.sub, aes(x = XXXX, y = YYYY, label = "+"), col = "white", size = max(1, floor(text.size/2)))
  }

  p

}




#' Description: Draw regression curve with smoothed error bars 
#' based on the Visually-Weighted Regression by Solomon M. Hsiang; see
#' http://www.fight-entropy.com/2012/07/visually-weighted-regression.html
#' The R implementation is based on Felix Schonbrodt's code (under MIT/FreeBSD license) from 
#' http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/
#'
#' Arguments:
#' @param formula formula
#' @param data data
#' @param title title
#' @param B number bootstrapped smoothers
#' @param shade plot the shaded confidence region?
#' @param shade.alpha should the CI shading fade out at the edges? (by reducing alpha; 0 = no alpha decrease, 0.1 = medium alpha decrease, 0.5 = strong alpha decrease)
#' @param spag plot spaghetti lines?
#' @param mweight should the median smoother be visually weighted?
#' @param show.lm should the linear regresison line be plotted?
#' @param show.median show median smoother
#' @param median.col median color
#' @param show.CI should the 95\% CI limits be plotted?
#' @param method the fitting function for the spaghettis; default: loess
#' @param bw define a default b/w-palette (TRUE/FALSE)
#' @param slices number of slices in x and y direction for the shaded region. Higher numbers make a smoother plot, but takes longer to draw. I wouldn'T go beyond 500
#' @param palette provide a custom color palette for the watercolors
#' @param ylim restrict range of the watercoloring
#' @param quantize either "continuous", or "SD". In the latter case, we get three color regions for 1, 2, and 3 SD (an idea of John Mashey)
#' @param ... further parameters passed to the fitting function, in the case of loess, for example, "span = .9", or "family = 'symmetric'"
#' @param verbose print information during execution
#'
#' Returns:
#' @return ggplot2 object
#'
#' @examples N <- 10; df <- data.frame(age = sort(runif(N, 0, 100)), hitchip = rnorm(N)); p <- vwReg(hitchip~age, df, shade = TRUE, mweight = TRUE, verbose = FALSE)
#'
#' @import ggplot2 plyr reshape2
#'
#' @export
#' @references See citation("microbiome") 
#' @author Based on the original version from Felix Schonbrodt. Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

vwReg <- function(formula, data, title="", B=1000, shade=TRUE, shade.alpha=.1, spag=FALSE, mweight=TRUE, show.lm=FALSE, show.median = TRUE, median.col = "white", show.CI=FALSE, method=loess, bw=FALSE, slices=200, palette=colorRampPalette(c("#FFEDA0", "#DD0000"), bias=2)(20), ylim=NULL, quantize = "continuous",  verbose = FALSE, ...) {

  # Circumvent global variable binding warnings
  x <- NA
  y <- NA
  dens.scaled <- NA
  alpha.factor <- NA
  value <- NA
  group <- NA
  M <- NA
  w3 <- NA
  UL <- NA
  LL <- NA


  IV <- all.vars(formula)[2]
  DV <- all.vars(formula)[1]
  data <- na.omit(data[order(data[, IV]), c(IV, DV)]) 
  if (bw) palette <- colorRampPalette(c("#EEEEEE", "#999999", "#333333"), bias=2)(20)
  if (verbose) {message("Computing boostrapped smoothers ...")}

  newx <- data.frame(seq(min(data[, IV]), max(data[, IV]), length=slices))
  colnames(newx) <- IV
  l0.boot <- matrix(NA, nrow=nrow(newx), ncol=B)

  l0 <- method(formula, data)

  for (i in 1:B) {
    data2 <- data[sample(nrow(data), replace=TRUE), ]
    data2 <- data2[order(data2[, IV]), ]

    if (class(l0)=="loess") {
      m1 <- method(formula, data2, control = loess.control(surface = "i", statistics="a", trace.hat="a"), ...)
    } else {
      m1 <- method(formula, data2, ...)
    }
      l0.boot[, i] <- predict(m1, newdata=newx)
    }

    # compute median and CI limits of bootstrap
    InstallMarginal("plyr")
    InstallMarginal("reshape2")
    CI.boot <- adply(l0.boot, 1, function(x) quantile(x, prob=c(.025, .5, .975, pnorm(c(-3, -2, -1, 0, 1, 2, 3))), na.rm=TRUE))[, -1]
    colnames(CI.boot)[1:10] <- c("LL", "M", "UL", paste0("SD", 1:7))
    CI.boot$x <- newx[, 1]
    CI.boot$width <- CI.boot$UL - CI.boot$LL

    # scale the CI width to the range 0 to 1 and flip it (bigger numbers = narrower CI)
    CI.boot$w2 <- (CI.boot$width - min(CI.boot$width))
    CI.boot$w3 <- 1-(CI.boot$w2/max(CI.boot$w2))

    # convert bootstrapped spaghettis to long format
    b2 <- melt(l0.boot)
    b2$x <- newx[,1]
    colnames(b2) <- c("index", "B", "value", "x")

    InstallMarginal("ggplot2")
    InstallMarginal("RColorBrewer")

    p1 <- ggplot(data, aes_string(x=IV, y=DV)) + theme_bw()

    if (shade) {
      quantize <- match.arg(quantize, c("continuous", "SD"))

        if (quantize == "continuous") {

	  if (verbose) {message("Computing density estimates for each vertical cut ...")}
	  flush.console()
  	    if (is.null(ylim)) {
	      min_value <- min(min(l0.boot, na.rm=TRUE), min(data[, DV], na.rm=TRUE))
	      max_value <- max(max(l0.boot, na.rm=TRUE), max(data[, DV], na.rm=TRUE))
	      ylim <- c(min_value, max_value)
	    }

	    # vertical cross-sectional density estimate
	    d2 <- ddply(b2[, c("x", "value")], .(x), function(df) {
	      res <- data.frame(density(df$value, na.rm=TRUE, n=slices, from=ylim[1], to=ylim[2])[c("x", "y")])

	      colnames(res) <- c("y", "dens")
	      return(res)
	    })
	    #}, .progress="text")

	    maxdens <- max(d2$dens)
	    mindens <- min(d2$dens)
	    d2$dens.scaled <- (d2$dens - mindens)/maxdens

	    ## Tile approach
	    d2$alpha.factor <- d2$dens.scaled^shade.alpha
	    p1 <- p1 + geom_tile(data=d2, aes(x=x, y=y, fill=dens.scaled, alpha=alpha.factor)) + scale_fill_gradientn("dens.scaled", colours=palette) + scale_alpha_continuous(range=c(0.001, 1))
	    }

	    if (quantize == "SD") {

	      ## Polygon approach
	      SDs <- melt(CI.boot[, c("x", paste0("SD", 1:7))], id.vars="x")
	      count <- 0
	      d3 <- data.frame()
	      col <- c(1,2,3,3,2,1)

	      for (i in 1:6) {
	        seg1 <- SDs[SDs$variable == paste0("SD", i), ]
		seg2 <- SDs[SDs$variable == paste0("SD", i+1), ]
		seg <- rbind(seg1, seg2[nrow(seg2):1, ])
		seg$group <- count
		seg$col <- col[i]
		count <- count + 1
		d3 <- rbind(d3, seg)
 	     }

	     p1 <- p1 + geom_polygon(data=d3, aes(x=x, y=value, color=NULL, fill=col, group=group)) + scale_fill_gradientn("dens.scaled", colours=palette, values=seq(-1, 3, 1))
	  }
	}	  

	if (verbose) {message("Build ggplot figure ...")}
	flush.console()
	if (spag) {
	  p1 <- p1 + geom_path(data=b2, aes(x=x, y=value, group=B), size=0.7, alpha=10/B, color="darkblue")
	}

	if (show.median) {
	  if (mweight) {
	    p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=M, alpha=w3^3), size=.6, linejoin="mitre", color=median.col)
	  } else {
	    p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=M), size = 0.6, linejoin="mitre", color=median.col)
	  }
	}

	# Confidence limits
	if (show.CI) {
  	  p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=UL, group=B), size=1, color="red")
  	  p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=LL, group=B), size=1, color="red")
	}

  # plain linear regression line
  if (show.lm) {p1 <- p1 + geom_smooth(method="lm", color="darkgreen", se=FALSE)}

  p1 <- p1 + geom_point(size=1, shape=21, fill="white", color="black")

  if (title != "") {
    p1 <- p1 + ggtitle(title)
  }

  p1 

}




#' Description: Project high-dimensional data on two-dimensional plane by various methods
#' 
#' Arguments:
#'   @param amat data matrix (samples x features)
#'   @param type projection type (options: PCA, MDS.classical, MDS.nonmetric, Sammon)
#' Returns:
#'   @return projected data matrix
#'
#' @export
#' @import MASS
#'
#' @examples data(peerj32); xy <- project.data(peerj32$microbes[,1:3])
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

project.data <- function (amat, type = "PCA") {

  if (type == "PCA") {
    if (nrow(amat) < ncol(amat)) {

      message("More samples than features, using sparse PCA")

      ## Spca example: we are selecting 50 variables on each of the PCs
      InstallMarginal("mixOmics")

      result <- mixOmics::spca(amat, ncomp = 2, center = TRUE, scale = TRUE, keepX = rep(50, 2))
      scores <- result$x
    } else {
      message("PCA")
      pca <- princomp(amat) # Classical PCA
      scores <- pca$scores
    }
    tab <- data.frame(scores[,1:2])
    rownames(tab) <- rownames(amat)
  } else if (type == "Sammon") {

    InstallMarginal("MASS")

    d <- as.dist(1-cor(t(amat)))
    # This gave the clearest visualization. 
    # Tuning magic parameter could still improve. 
    # Try for instance magic = 0.05.
    fit <- sammon(d, k = 2) 
    # Plot solution 
    tab <- data.frame(list(Comp.1 = fit$points[,1], Comp.2 = fit$points[,2]))
    rownames(tab) <- rownames(amat)
  } else if (type == "MDS.classical") {
    d <- as.dist(1-cor(t(amat)))
    fit <- cmdscale(d, eig=TRUE, k=2) # classical MDS
    tab <- data.frame(list(Comp.1 = fit$points[,1], Comp.2 = fit$points[,2]))    
  } else if (type == "MDS.nonmetric") {
    d <- as.dist(1-cor(t(amat)))
    fit <- isoMDS(d, k=2)             # nonmetric MDS
    tab <- data.frame(list(Comp.1 = fit$points[,1], Comp.2 = fit$points[,2]))    
  }  

  # TODO Kernel-PCA
  #library(kernlab)
  #kpc <- kpca(~., data=as.data.frame(x.train), kernel="rbfdot", features = 2)
  #Print the principal component vectors
  #pcv(kpc)
  #Plot the data projection on the components
  #par(mfrow=c(2,2))
  #plot(rotated(kpc), col = as.integer(as.factor(ann[rownames(x.train),"time"])), xlab="1st Principal Component", ylab="2nd Principal Comp  onent")
  #plot(rotated(kpc), col = as.integer(as.factor(ann[rownames(x.train),"lipids.group"])), xlab="1st Principal Component", ylab="2nd Principal Component")
  #embed remaining points 
  #emb <- predict(kpc, x.test)
  #plot(rotated(kpc), col = as.integer(as.factor(ann[rownames(x.train),"lipids.group"])), xlab="1st Principal Component", ylab="2nd Principal Component")
  #points(emb, col = as.integer(as.factor(ann[rownames(x.train),"lipids.group"])))
    
  colnames(tab) <- c("Comp.1","Comp.2")

  tab
}


#' Visualize a matrix with one or two-way color scale. 
#' TODO: one-way color scale
#'
#' Fast investigation of matrix objects; standard visualization choices are 
#' made automatically; fast and easy-to-use but does not necessarily provide 
#' optimal visualization.
#'
#' @param mat matrix
#' @param type String. Specifies visualization type. Options: "oneway" (color scale ranges from white to dark red; the color can be changed if needed); "twoway" (color scale ranges from dark blue through white to dark red; colors can be changed if needed)
#' @param midpoint middle point for the color plot: smaller values are shown with blue, larger are shown with red in type = "twoway"
#' @param palette Optional. Color palette.
#' @param colors Optional. Colors.
#' @param col.breaks breakpoints for the color palette
#' @param interval interval for palette color switches
#' @param plot.axes String. Indicates whether to plot x-axis ("x"), y-axis ("y"), or both ("both").
#' @param row.tick interval for plotting row axis texts
#' @param col.tick interval for plotting column axis texts
#' @param cex.xlab use this to specify distinct font size for the x axis
#' @param cex.ylab use this to specify distinct font size for the y axis
#' @param xlab optional x axis labels
#' @param ylab optional y axis labels
#' @param limit.trunc color scale limit breakpoint
#' @param mar image margins
#' @param ... optional parameters to be passed to function 'image', see help(image) for further details
#' @return A list with the color palette (colors), color breakpoints (breaks), and palette function (palette.function)
#' @export
#' 
#' @references See citation("microbiome") 
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples mat <- rbind(c(1,2,3,4,5), c(1, 3, 1), c(4,2,2)); PlotMatrix(mat, "twoway", midpoint = 3) 
#' @keywords utilities

PlotMatrix <- function (mat, type = "twoway", midpoint = 0, 
	      	        palette = NULL, colors = NULL, col.breaks = NULL, interval = .1, 
			plot.axes = "both",
			row.tick = 1, col.tick = 1, 
			cex.xlab = .9, cex.ylab = .9, 
			xlab = NULL, ylab = NULL,
			limit.trunc = 0, mar = c(5, 4, 4, 2), ...) {

  # Center the data and color breakpoints around the specified midpoint
  mat <- mat - midpoint

  if (length(col.breaks) == 0)  {
    m <- max(round(max(abs(mat)), limit.trunc) - interval, 0)

    mm <- m + interval/2

    vals <- seq(interval/2, mm, interval)

    # Set col.breaks evenly around zero
    col.breaks  <- c(-(m + 1e6), c(-rev(vals), vals), m+1e6)
  }

  if (is.null(palette)) {
    my.palette <- colorRampPalette(c("blue", "white", "red"), space = "rgb")
  } else if (palette == "blue-black-red") {
    my.palette <- colorRampPalette(c("blue", "black", "red"), space = "rgb")
  } else if (palette == "blue-white-red") {
    my.palette <- colorRampPalette(c("blue", "white", "red"), space = "rgb")
  } else if (palette == "blue-white-yellow") {
    my.palette <- colorRampPalette(c("blue", "white", "yellow"), space = "rgb")
  } else if (palette == "blue-black-yellow") {
    my.palette <- colorRampPalette(c("blue", "black", "yellow"), space = "rgb")
  } else if (palette == "bw") {
    gray.palette <- function (int) {
      gray(seq(0,1,length=int))
    }
    my.palette <- gray.palette
  }
  
  # if mycolors is provided it overrides palette
  if (is.null(colors)) { colors <- my.palette(length(col.breaks) - 1) }

  # transpose and revert row order to plot matrix in the same way it
  # appears in its numeric form
  par(mar = mar)
  
  matm <- matrix(mat[rev(seq(nrow(mat))),], ncol = ncol(mat))
  dimnames(matm) <- dimnames(mat)
  mat <- matm

  nsamples <- ncol(mat)
  nfeats <- nrow(mat)
  if (nfeats == 1) {mat <- mat}
  image(t(mat), col = colors, xaxt = 'n', yaxt = 'n', zlim = range(col.breaks), breaks = col.breaks, ...)

  if (plot.axes == "both" || plot.axes == TRUE) {

    if (is.null(xlab)) {
      v <- seq(1, nsamples, col.tick) # take every nth index
      axis(1, at = seq(0,1,length = nsamples)[v], labels = colnames(mat)[v], cex.axis=cex.xlab, las=2, ...)    
    } else {
      axis(1, at = seq(0,1,length = nsamples), labels = xlab, cex.axis=cex.xlab, las=2, ...)    
    }

    if (is.null(ylab)) {
      v <- seq(1, nfeats, row.tick) # take every nth index
      axis(2, at = seq(0,1,length = nfeats)[v], labels = rev(rownames(mat))[v], cex.axis=cex.ylab, las=2, ...)
    } else {  
      axis(2, at = seq(0,1,length = nfeats), labels = ylab, cex.axis=cex.ylab, las=2, ...)
    }

  } else if (plot.axes == "x") {

    if (is.null(xlab)) {
      v <- seq(1, nsamples, col.tick) # take every nth index
      axis(1, at = seq(0,1,length = nsamples)[v], labels = colnames(mat)[v], cex.axis=cex.ylab, las=2)    
    } else {
      axis(1, at = seq(0,1,length = nsamples), labels = ylab, cex.axis=cex.ylab, las=2)    
    }

  } else if (plot.axes == "y") {

    if (is.null(ylab)) {
      v <- seq(1, nfeats, row.tick) # take every nth index
      axis(2, at = seq(0, 1, length = nfeats)[v], labels = rev(rownames(mat))[v], cex.axis = cex.xlab, las = 2)
    } else {  
      axis(2, at = seq(0, 1, length = nfeats), labels = ylab, cex.axis=cex.xlab, las=2)
    }
  }
  
  # Return default margins
  par(mar = c(5, 4, 4, 2) + 0.1)
 
  return(list(colors = colors, breaks = col.breaks + midpoint, palette.function = my.palette))
      	  
}






#' htree.plot
#' Description: Plot hierarchical clustering for the input data in absolute
#' and log10 scale using euclidean and pearson correlation similarities. 
#' Intended for internal use in the run.profiling.script function. 
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param method hierarchical clustering method
#'   @param metric clustering similarity measure
#' Returns:
#'   @return NULL
#'
#' @export
#' @examples data(peerj32); tmp <- htree.plot(peerj32$microbes[,1:5])
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

htree.plot <- function (dat, method = "complete", metric = "pearson") {

  # Plot CLUSTER TREES TO A GRAPHICS WINDOW
  # Euclidean & Correlation / Raw & Log10

  if(ncol(dat) > 2){

    if (metric == "euclidean") {
      # Metric: Euclidean
      hc <- hclust(dist(t(dat)), method = method)
      #hc.raw.eu <- hclust(dist(t(dat)), method = method)
      #hc.log10.eu <- hclust(dist(t(log10(dat + 1))), method = method)
    } else if (metric %in% c("pearson", "spearman")) {

      # 'Metric': Correlation
      hc <- hclust(as.dist(1 - cor(dat, use = "pairwise.complete.obs", method = metric)), method = method)

      #hc.raw.cor <- hclust(as.dist(1 - cor(dat, use = "complete.obs")), method = method)
      #hc.log10.cor <- hclust(as.dist(1 - cor(log10(dat + 1), use = "complete.obs")), method = method)
    }


    # Plot all hclust trees in a single figure
    #x11() 
    #par(mfrow=c(2,2))
    #plot(hc.log10.eu, hang = -1, main = "hclust/euclid/oligo/log10", xlab = "Samples")
    #plot(hc.raw.eu, hang = -1, main = "hclust/euclid/oligo/raw", xlab = "Samples")
    #plot(hc.log10.cor, hang = -1, main = "hclust/pearson/oligo/log10", xlab = "Samples")
    #plot(hc.raw.cor, hang = -1, main = "hclust/pearson/oligo/raw", xlab = "Samples")
    plot(hc, hang = -1, main = paste("hclust/", metric, sep = ""), xlab = "Samples")

  } else {
    warning("Three or more samples required for clustering - skipped.\n")
  }

  return(NULL)

}




#' phylo.barplot
#'
#' Description: Barplot for *ITChip sample (across taxa) with higher-level taxonomic groups indicated by colours.
#'
#' Arguments:
#'   @param x Data vector across taxa (each element should be named by taxon)
#'   @param color.level Higher-order phylogenetic level to indicate by colors
#'   @param phylogeny.info oligo-phylotype mappings
#'   @param title title
#'   @param plot draw plot TRUE/FALSE
#'   @param sort sort the effects by magnitude
#'
#' Returns:
#'   @return ggplot2 object
#'
#' @export
#' @examples # NOT RUN: phylogeny.info <- GetPhylogeny("HITChip", "filtered"); signal <-  unlist(peerj32$microbes[1, 1:10]); p <- phylo.barplot(signal, color.level = "L1", phylogeny.info = phylogeny.info)
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

phylo.barplot <- function (x, color.level = "L1", phylogeny.info = NULL, title = NULL, plot = TRUE, sort = TRUE) {

  if (is.null(phylogeny.info)) {
    warning("phylogeny.info not specified, assuming HITChip phylogeny")
    phylogeny.info <- GetPhylogeny("HITChip", "filtered")
  }

  taxa <- names(x)

  for (tax.lev in c("oligoID", "species", "L1", "L2")) {
    if (all(taxa %in% phylogeny.info[[tax.lev]])) {x.level <- tax.lev}
  }

  # Collect all into a data.frame
  df <- data.frame(list(taxa = taxa))
  
  # Assign higher-level taxonomic groups
  df[[color.level]] <- unlist(droplevels(levelmap(taxa, level.from = x.level, level.to = color.level, phylogeny.info = phylogeny.info)))
  df[["color.level"]] <- df[[color.level]]

  df[["x"]] <- x

  # Define colors for L1/L2 groups
  colors <- rainbow(length(unique(df[[color.level]])))
  names(colors) <- as.character(unique(df[[color.level]]))

  # Rearrange data.frame
  m <- melt(df)

  # Sort by x (ie. change order of factors for plot)
  if (sort) {
    df <- within(df, taxa <- factor(taxa, levels = taxa[order(abs(x))]))
  }

  # Plot the image
  p <- ggplot(aes(x = taxa, y = x, fill = color.level), data = df)
  p <- p + scale_fill_manual(values = colors[as.character(levels(df[[color.level]]))])
  p <- p + geom_bar(position = "identity", stat = "identity") + theme_bw() + coord_flip()
  p <- p + ylab("Signal") + xlab("") + ggtitle(title)
  p <- p + theme(legend.position = "right")
  p <- p + theme(panel.border = element_rect())

  if (plot) {
    print(p)
  }

  p

}
