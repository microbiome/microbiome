# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Draw regression curve with smoothed error bars 
#' based on the Visually-Weighted Regression by Solomon M. Hsiang; see
#' http://www.fight-entropy.com/2012/07/visually-weighted-regression.html
#' The R implementation is based on Felix Schönbrodt's code from 
#' http://www.nicebread.de/visually-weighted-regression-in-r-a-la-solomon-hsiang/
#'
#' Arguments:
#'   @param formula formula
#'   @param data data
#'   @param B number bootstrapped smoothers
#'   @param shade: plot the shaded confidence region?
#'   @param spag: plot spaghetti lines?
#'   @param mweight: should the median smoother be visually weighted?
#'   @param show.lm: should the linear regresison line be plotted?
#'   @param show.CI: should the 95% CI limits be plotted?
#'
#' Returns:
#'   @return ggplot2 object
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

vwReg <- function(formula, data, B=1000, shade=TRUE, spag=FALSE, mweight=TRUE, show.lm=FALSE, show.CI=FALSE) {

#    Compute smoothers from 1000 bootstrap samples of the original sample (this results in a spaghetti plot)
#    Calculate a density estimate for each vertical cut through the bootstrapped smoothers. The area under the density curve always is 1, so the ink is constant for each y-slice.
#    Shade the figure according to these density estimates.
# -> Implemented in p4 below
## build a demo data set
#set.seed(1)
#x <- rnorm(200, 0.8, 1.2) 
#e <- rnorm(200, 0, 2)*(abs(x)^1.5 + .5)
#y <- 8*x - x^3 + e
#df <- data.frame(x, y)
 
#p1 <- vwReg(y~x, df)
#p2 <- vwReg(y~x, df, shade=FALSE, spag=TRUE)
#p3 <- vwReg(y~x, df, shade=TRUE, spag=FALSE, mweight=TRUE, show.CI=TRUE, show.lm=TRUE)
#p4 <- vwReg(y~x, df, shade=FALSE, spag=TRUE, show.lm=TRUE)

	IV <- all.vars(formula)[2]
	DV <- all.vars(formula)[1]
	data <- na.omit(data[order(data[, IV]), c(IV, DV)])
 
	print("Computing boostrapped smoothers ...")
	steps <- nrow(data)*3	# three times more prediction points on the x-axis than original data points - for smoother smoothers
	newx <- seq(min(data[, IV]), max(data[, IV]), length=steps)
	l0.boot <- matrix(NA, nrow=steps, ncol=B)
	for (i in 1:B) {
		data2 <- data[sample(nrow(data), replace=TRUE), ]
		data2 <- data2[order(data2[, IV]), ]
		l0.boot[, i] <- predict(loess(formula, data2), newdata=newx)
	}
 
	# compute median and CI limits of bootstrap
	library(plyr)
        library(reshape2)
	CI.boot <- adply(l0.boot, 1, function(x) quantile(x, prob=c(.025, .5, .975), na.rm=TRUE))[, -1]
	colnames(CI.boot)[1:3] <- c("LL", "M", "UL")
	CI.boot$x <- newx
	CI.boot$width <- CI.boot$UL - CI.boot$LL
 
	# scale the CI width to the range 0 to 1 and flip it (bigger numbers = narrower CI)
	CI.boot$w2 <- (CI.boot$width - min(CI.boot$width))
	CI.boot$w3 <- 1-(CI.boot$w2/max(CI.boot$w2))
 
	# convert bootstrapped spaghettis to long format
	b2 <- melt(l0.boot)
	b2$x <- newx
	colnames(b2) <- c("index", "B", "value", "x")
 
	print("ggplot prints the figure ...")
	library(ggplot2)
	library(RColorBrewer)
 
	p1 <- ggplot(data, aes_string(x=IV, y=DV)) + theme_bw()
 
 
	if (shade == TRUE) {
		print("Computing density estimates for each vertical cut ...")
		# vertical cross-sectional density estimate
		d2 <- ddply(b2[, c("x", "value")], .(x), function(df) {
			res <- data.frame(density(df$value, na.rm=TRUE, n=100)[c("x", "y")])
			colnames(res) <- c("y", "dens")
			return(res)
		}, .progress="text")
 
		maxdens <- max(d2$dens)
		mindens <- min(d2$dens)
		d2$dens.scaled <- (d2$dens - mindens)/maxdens	
 
		# alpha scaling
		#p1 <- p1 + geom_point(data=d2, aes(x=x, y=y, alpha=dens.scaled), size=0.4, color="red")
 
		# color scaling
		p1 <- p1 + geom_point(data=d2, aes(x=x, y=y, color=dens.scaled), size=1.4) + scale_color_gradientn("dens.scaled", colours=brewer.pal(9, "YlOrRd"))
	}
 
	if (spag==TRUE) {
		library(reshape2)
		p1 <- p1 + geom_path(data=b2, aes(x=x, y=value, group=B), size=0.7, alpha=15/B, color="darkblue")
	}
 
	if (mweight == TRUE) {
		p1 <- p1 + geom_point(data=CI.boot, aes(x=x, y=M, alpha=w3), size=1, linejoin="mitre", color="darkred")
	} else {
		p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=M), size = 0.3, linejoin="mitre", color="darkred")
	}
 
	# Confidence limits
	if (show.CI == TRUE) {
		#p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=UL, group=B), size=1, color="red")
		#p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=LL, group=B), size=1, color="red")

		p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=UL), size=1, color="red")
		p1 <- p1 + geom_path(data=CI.boot, aes(x=x, y=LL), size=1, color="red")
	}
 
	p1 <- p1 + geom_point(size=1)
 
	# plain linear regression line
	if (show.lm==TRUE) {p1 <- p1 + geom_smooth(method="lm", color="darkgreen", se=FALSE)}
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
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

project.data <- function (amat, type = "PCA") {

  if (type == "PCA") {
    if (nrow(amat) < ncol(amat)) {

      message("More samples than features, using sparse PCA")
      ## Spca example: we are selecting 50 variables on each of the PCs
      library(mixOmics)
      result <- spca(amat, ncomp = 2, center = TRUE, scale. = TRUE, keepX = rep(50, 2))
      scores <- result$x
    } else {
      message("PCA")
      pca <- princomp(amat) # Classical PCA
      scores <- pca$scores
    }
    tab <- data.frame(scores[,1:2])
    rownames(tab) <- rownames(amat)
  } else if (type == "Sammon") {
    library(MASS)
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
#' @references See citation("microbiome") 
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # mat <- rbind(c(1,2,3,4,5), c(1, 3, 1), c(4,2,2)); PlotMatrix(mat, "twoway", midpoint = 3) 
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
  image(t(mat[rev(seq(nrow(mat))),]), col = colors, xaxt = 'n', yaxt = 'n', zlim = range(col.breaks), breaks = col.breaks, ...)

  if (plot.axes == "both" || plot.axes == TRUE) {

    if (is.null(xlab)) {
      v <- seq(1, ncol(mat), col.tick) # take every nth index
      axis(1, at = seq(0,1,length = ncol(mat))[v], labels = colnames(mat)[v], cex.axis=cex.xlab, las=2, ...)    
    } else {
      axis(1, at = seq(0,1,length = ncol(mat)), labels = xlab, cex.axis=cex.xlab, las=2, ...)    
    }

    if (is.null(ylab)) {
      v <- seq(1, nrow(mat), row.tick) # take every nth index
      axis(2, at = seq(0,1,length = nrow(mat))[v], labels = rev(rownames(mat))[v], cex.axis=cex.ylab, las=2, ...)
    } else {  
      axis(2, at = seq(0,1,length = nrow(mat)), labels = ylab, cex.axis=cex.ylab, las=2, ...)
    }

  } else if (plot.axes == "x") {

    if (is.null(xlab)) {
      v <- seq(1, ncol(mat), col.tick) # take every nth index
      axis(1, at = seq(0,1,length = ncol(mat))[v], labels = colnames(mat)[v], cex.axis=cex.ylab, las=2)    
    } else {
      axis(1, at = seq(0,1,length = ncol(mat)), labels = ylab, cex.axis=cex.ylab, las=2)    
    }

  } else if (plot.axes == "y") {

    if (is.null(ylab)) {
      v <- seq(1, nrow(mat), row.tick) # take every nth index
      axis(2, at = seq(0, 1, length = nrow(mat))[v], labels = rev(rownames(mat))[v], cex.axis = cex.xlab, las = 2)
    } else {  
      axis(2, at = seq(0, 1, length = nrow(mat)), labels = ylab, cex.axis=cex.xlab, las=2)
    }
  }
  
  # Return default margins
  par(mar = c(5, 4, 4, 2) + 0.1)
 
  return(list(colors = colors, breaks = col.breaks + midpoint, palette.function = my.palette))
      	  
}




#' Description: ggplot2::ggplot2 theme
#'              
#' Arguments:
#'   @param colour colour
#'   @param size size
#'   @param linetype linetype
#' Returns:
#'   @return set theme
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

theme_bottom_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot2::ggplot(...) + opts( panel.border=theme_bottom_border() ) + ...
  structure(
    function(x = 0, y = 0, width = 1, height = 1, ...) {
      polylineGrob(
        x=c(x, x+width), y=c(y,y), ..., default.units = "npc",
        gp=gpar(lwd=size, col=colour, lty=linetype),
      )
    },
    class = "theme",
    type = "box",
    call = match.call()
  )
}



#' plot.htrees
#' Description: Plot hierarchical clustering for the input data in absolute
#' and log10 scale using euclidean and pearson correlation similarities. 
#' Intended for internal use in the run.profiling.script function. 
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param method hierarchical clustering method
#' Returns:
#'   @return NULL
#'
#' @export
#' @examples # dat <- read.profiling(params$wdir, "species", "rpa"); tmp <- plot.htrees(dat)
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plot.htrees <- function (dat, method = "complete") {

  # Plot CLUSTER TREES TO A GRAPHICS WINDOW
  # Euclidean & Correlation / Raw & Log10

  if(ncol(dat) > 2){

    # Metric: Euclidean
    hc.raw.eu <- hclust(dist(t(dat)), method = method)
    hc.log10.eu <- hclust(dist(t(log10(dat + 1))), method = method)

    # 'Metric': Correlation
    hc.raw.cor <- hclust(as.dist(1 - cor(dat, use = "complete.obs")), method = method)
    hc.log10.cor <- hclust(as.dist(1 - cor(log10(dat + 1), use = "complete.obs")), method = method)

    # Plot all hclust trees in a single figure
    x11() 
    par(mfrow=c(2,2))
    plot(hc.log10.eu, hang = -1, main = "hclust/euclid/oligo/log10", xlab = "Samples")
    plot(hc.raw.eu, hang = -1, main = "hclust/euclid/oligo/raw", xlab = "Samples")
    plot(hc.log10.cor, hang = -1, main = "hclust/pearson/oligo/log10", xlab = "Samples")
    plot(hc.raw.cor, hang = -1, main = "hclust/pearson/oligo/raw", xlab = "Samples")

  } else {
    warning("Three or more samples required for clustering - skipped.\n")
  }

  return(NULL)

}


#' add.heatmap
#' Description: Add oligprofile heatmap into output directory
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param output.dir output data directory
#'   @param output.dir output file name
#'   @param phylogeny.info oligo-phylotype mappings
#'   @param ppcm figure size
#'   @param hclust.method hierarchical clustering method
#'   @param palette color palette ("white/black" / "white/blue" / "black/yellow/white")
#'   @param level taxonomic level to show
#'   @param metric clustering metric
#'   @param figureratio figure ratio
#'   @param fontsize font size
#'   @param tree.display tree.display
#'
#' Returns:
#'   @return Plotting parameters
#'
#' @export
#' @examples # dat <- read.profiling(params$wdir, "species", "rpa"); hc <- add.heatmap(dat)
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

add.heatmap <- function (dat, output.dir, output.file = NULL, phylogeny.info, ppcm = 150, 
	         hclust.method = "complete", palette = "white/black", level = "L1", metric = "pearson", 
  		 figureratio = 10, fontsize = 40, tree.display = TRUE) {

  # dat <- finaldata[["oligo"]]; output.dir = params$wdir;  output.file = NULL; phylogeny.info = phylogeny.info; ppcm = 150; hclust.method = "complete"; palette = "white/blue"; level = "L2"; metric = "pearson"; figureratio = 12; fontsize = 12; tree.display = TRUE

  if (is.null(output.file)) {
    output.file <- paste(output.dir,"/", gsub(" ", "", level), "-oligoprofileClustering.png",sep="")
  }		 

  hc.params <- list()
  if( ncol(dat) >= 3 ) {

    message(paste("Storing oligo heatmap in", output.file))  
    hc.params$ppcm <- ppcm
    hc.params$output.file <- output.file

    # PLOT THE HEATMAP
    # figure width as a function of the number of the samples
    plotdev <- png(filename = output.file, 
  	    width = max(trunc(ppcm*21), trunc(ppcm*21*ncol(dat)/70)), 
	    height = trunc(ppcm*29.7)) 
    try(hc.params <- PlotPhylochipHeatmap(data = dat,
                phylogeny.info = phylogeny.info,
                metric = metric,
                level = level,
                tree.display = tree.display,
                palette = palette,
                fontsize = fontsize,
                figureratio = figureratio, 
		hclust.method = hclust.method)) 

    dev.off()
  }

  hc.params

}


#' PlotPhylochipHeatmap
#'
#' Description: Plots heatmap of the oligo profiles together with phylotype grouping and sample clustering
#'
#' Arguments:
#'   @param data oligoprofile data in original (non-log) domain
#'   @param phylogeny.info oligo-phylotype mappings
#'   @param metric clustering metric
#'   @param level taxonomic level to show (L0 / L1 / L2 / species)
#'   @param tree.display tree.display
#'   @param palette color palette ("white/black" / "white/blue" / "black/yellow/white")
#'   @param fontsize font size
#'   @param figureratio figure ratio
#'   @param hclust.method hierarchical clustering method
#'
#' Returns:
#'   @return parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PlotPhylochipHeatmap <- function (data,
                         phylogeny.info,
                         metric = "pearson", 
                         level = "L1", 
                         tree.display = TRUE, 
                         palette = "white/black", 
                         fontsize = 40, 
                         figureratio = 10, 
			 hclust.method = "complete") {

  # data = dat; metric = "pearson"; level = "L2"; tree.display = TRUE; palette = "white/black"; fontsize = 40; figureratio = 10; hclust.method = "complete"

  params <- c(metric = metric, level = level, tree.display = tree.display, palette = palette, 
  	      fontsize = fontsize, figureratio = figureratio, hclust.method = hclust.method)
			 
  if (is.character(palette)) {
    palette  <- list.color.scales()[[palette]]
  }
              
   par(ps = fontsize, xpd = NA)
   paper <- par("din")

   if (level == "oligo") { level <- "oligoID" }
   tax.order <- order(as.character(phylogeny.info[[level]]), na.last = FALSE)

   nainds <- is.na(phylogeny.info[, level])
   if (sum(nainds) > 0) {
     phylogeny.info[nainds, level] <- '-'  # replace empty maps
   }

   levs <- unlist(lapply(split(phylogeny.info[[level]], as.factor(phylogeny.info[[level]])), length))
   # order the rows in phylogeny.info by level
   phylogeny.info <- phylogeny.info[tax.order,]
   phylogeny.info <- phylogeny.info[phylogeny.info$oligoID %in% rownames(data), ]

   annwidth <- max(strwidth(names(levs),units="inch"))*2.54*1.2
   profilewidth <- max(strheight(names(levs),units="inch"))*2.54*dim(data)[2]*1.6
   figureheight <- paper[2]*2.54*0.9

   # prevent outliers from determining the ends of the colour scale
   limits <- quantile(data, c(0.001,0.999), na.rm = TRUE)
   limits <- limits*c(0.98, 1.02)

   # calculate clustering based on oligoprofile
   if (metric == "euclidean") {
    hc <- hclust(dist(t(data)), method = hclust.method)
   } else if (metric == "pearson") {
    hc <- hclust(as.dist(1 - cor(data, use = "pairwise.complete.obs")), method = hclust.method)
   }

   data[data < limits[1]] <- limits[1]
   data[data > limits[2]] <- limits[2]
   if (!is.na(figureratio)) {
      heights = c(figureratio/100, (100-figureratio)/100)
   } else {
      if (tree.display) {
         heights = c(15/100, 85/100)
      } else {
         heights = c(6/100, 94/100)
      }
   }

   if (tree.display) {

      data <- data[,hc$order] 
      layout(matrix(c(3,0,1,2),ncol=2,byrow=TRUE),widths=lcm(c(profilewidth,annwidth)),heights=lcm(figureheight*heights))

   } else {    

      layout(matrix(c(3,0,1,2),ncol=2,byrow=TRUE),widths=lcm(c(profilewidth,annwidth)),heights=lcm(figureheight*heights))

   }

   par(mar = c(1,1,0,0), oma = c(0,0,0,0))
   
   img <- as.matrix(rev(as.data.frame(t(data[as.character(phylogeny.info$oligoID),]))))
   image(z=img, col=palette, axes=FALSE, frame.plot=TRUE, zlim=limits)
   plot.new()
   par(mar = c(1, 0, 0, 1), usr = c(0, 1, 0, 1), xaxs = 'i', yaxs = 'i')

   rect(xleft = rep(0,length(levs)),
        ybottom = c(0,cumsum(rev(levs))[1:length(levs)-1])/sum(levs),
        xright = rep(1,length(levs)),ytop=cumsum(rev(levs))/sum(levs), border = 'grey')

   text(x = c(0.03), y = (cumsum(rev(levs))-rev(levs/2))/sum(levs), labels = (names(rev(levs))), pos = 4)

   if (tree.display) {

      par(mar=c(0.2,1.5,1,0.5),usr=c(0,1,0,1))
      plot(hc, axes = FALSE, ann = FALSE, hang = -1)

   } else {

      plot.new()
      par(mar=c(0,1,1,0),usr=c(0,1,0,1),xaxs='i',yaxs='i')
      text(x=0.5/length(colnames(data))+(seq(along.with=colnames(data))-1)/length(colnames(data)),
      y = c(0.15), labels = colnames(data), pos = 3, cex = 0.8, srt = 90)

   }

  params

}


