# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' PlotPhylochipHeatmap
#'
#' Description: Plots heatmap of the oligo profiles together with phylotype grouping and sample clustering
#'
#' Arguments:
#'   @param data oligoprofile data in original (non-log) domain
#'   @param phylogeny oligo-phylotype mappings
#'   @param metric clustering metric
#'   @param tax.level taxonomic level to show
#'   @param include.tree include.tree
#'   @param palette color palette
#'   @param fontsize font size
#'   @param figureratio figure ratio
#'   @param hclust.method hierarchical clustering method
#'
#' Returns:
#'   @return plots the profile
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PlotPhylochipHeatmap <- function (data,
                         phylogeny,
                         metric = "euclidian", 
                         tax.level = "level 2", 
                         include.tree = FALSE, 
                         palette = "black/yellow/white",
                         fontsize = 12, 
                         figureratio = 12, 
			 hclust.method = "ward") {

  # metric = "euclidian"; tax.level = "level 2"; include.tree = FALSE; palette = "black/yellow/white"; fontsize = 12; figureratio = 12; hclust.method = "ward"
			 
  metric <- list.clustering.metrics()[[metric]]

  if (is.character(palette)) {
    palette  <- list.color.scales()[[palette]]
    # colorRampPalette(c("black","yellow","white"), bias=0.5, interpolate = "linear")(100)
   }
              
   par(ps = fontsize, xpd = NA)
   paper <- par("din")

   if (tax.level == "oligo") {tax.level <- "oligoID"}
   tax.order <- order(phylogeny[[tax.level]], na.last=FALSE)

   nainds <- is.na(phylogeny[,tax.level])
   if (sum(nainds) > 0) {
     phylogeny[nainds,tax.level] <- '-'  # replace empty maps
   }

   levs <- unlist(lapply(split(phylogeny[[tax.level]], as.factor(phylogeny[[tax.level]])), length))
   # order the rows in phylogeny by tax.level
   phylogeny <- phylogeny[tax.order,]

   annwidth = max(strwidth(names(levs),units="inch"))*2.54*1.2
   profilewidth = max(strheight(names(levs),units="inch"))*2.54*dim(data)[2]*1.6
   figureheight = paper[2]*2.54*0.9

   # prevent outliers from determining the ends of the colour scale
   limits = quantile(data, c(0.001,0.999),na.rm=TRUE)
   limits = limits*c(0.98,1.02)

   # calculate clustering based on oligoprofile
   if (metric == 'euclidian') {
    hc <- hclust(dist(t(data)), method=hclust.method)
   } else {
    hc <- hclust(as.dist(1-cor(data,use="pairwise.complete.obs")), method=hclust.method)
   }

   data[data<limits[1]]=limits[1]
   data[data>limits[2]]=limits[2]
   if (!is.na(figureratio)) {
      heights = c(figureratio/100, (100-figureratio)/100)
   } else {
      if (include.tree) {
         heights = c(15/100, 85/100)
      } else {
         heights = c(6/100, 94/100)
      }
   }

   if (include.tree) {

      data <- data[,hc$order] 
      layout(matrix(c(3,0,1,2),ncol=2,byrow=TRUE),widths=lcm(c(profilewidth,annwidth)),heights=lcm(figureheight*heights))

   } else {    

      layout(matrix(c(3,0,1,2),ncol=2,byrow=TRUE),widths=lcm(c(profilewidth,annwidth)),heights=lcm(figureheight*heights))

   }

   par(mar=c(1,1,0,0),oma=c(0,0,0,0))
   img <- as.matrix(rev(as.data.frame(t(data[phylogeny$oligoID,]))))
   image(z=img, col=palette, axes=FALSE, frame.plot=TRUE, zlim=limits)
   plot.new()
   par(mar=c(1,0,0,1),usr=c(0,1,0,1),xaxs='i',yaxs='i')

   rect(xleft = rep(0,length(levs)),
        ybottom = c(0,cumsum(rev(levs))[1:length(levs)-1])/sum(levs),
        xright = rep(1,length(levs)),ytop=cumsum(rev(levs))/sum(levs), border = 'grey')

   text(x = c(0.03), y = (cumsum(rev(levs))-rev(levs/2))/sum(levs), labels = (names(rev(levs))), pos = 4)

   if (include.tree) {

      par(mar=c(0.2,1.5,1,0.5),usr=c(0,1,0,1))
      plot(hc, axes = FALSE, ann = FALSE, hang = -1)

   } else {

      plot.new()
      par(mar=c(0,1,1,0),usr=c(0,1,0,1),xaxs='i',yaxs='i')
      text(x=0.5/length(colnames(data))+(seq(along.with=colnames(data))-1)/length(colnames(data)),
      y = c(0.15), labels = colnames(data), pos = 3, cex = 0.8, srt = 90)

   }

}

# -----------------------------------------------------------------------


#' Description: Project high-dimensional data on two-dimensional plane by various methods
#' 
#' Arguments:
#'   @param amat data matrix (samples x features)
#'   @param type projection type
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
    fit <- sammon(d, k = 2) # This gave the clearest visualization. Tuning magic parameter could still improve. Try for instance magic = 0.05.
    # Plot solution 
    tab <- data.frame(list(Comp.1 = fit$points[,1], Comp.2 = fit$points[,2]))
    rownames(tab) <- rownames(amat)
  } 

  # TODO MDS
  #fit <- cmdscale(d, eig=TRUE, k=2) # classical MDS
  #fit <- isoMDS(d, k=2)             # nonmetric MDS

#library(MASS)
#d <- as.dist(1-cor(t(amat)))
#fit <- cmdscale(d, eig=TRUE, k=2) # classical MDS
#fit <- isoMDS(d, k=2)             # nonmetric MDS
#fit <- sammon(d, k = 2) # This gave the clearest visualization. Tuning magic parameter could still improve. Try for instance magic = 0.05.
# Plot solution 
#x <- fit$points[,1]
#y <- fit$points[,2]
#s <- rownames(amat)
#s <- c(ns,rs)
#par(mar = c(5,5,1,1))
#plot(x[s], y[s], xlab="Coordinate 1", ylab="Coordinate 2", type="n", cex.lab = 1.8, las = 1, cex.axis = 1.3)
#lab <- gsub("FLN","",gsub("_","",annot[s,]$sampleID))
#char <- "W"; lab <- sapply(strsplit(lab, char), function (x) {if (length(x)>1) {gsub(" ","",paste(x[[1]], char, collapse = ""))} else {x}})
#char <- "R"; lab <- sapply(strsplit(lab, char), function (x) {if (length(x)>1) {gsub(" ","",paste(x[[1]], char, collapse = ""))} else {x}})
#char <- "S"; lab <- sapply(strsplit(lab, char), function (x) {if (length(x)>1) {gsub(" ","",paste(x[[1]], char, collapse = ""))} else {x}})
#char <- "M"; lab <- sapply(strsplit(lab, char), function (x) {if (length(x)>1) {gsub(" ","",paste(x[[1]], char, collapse = ""))} else {x}})
#char <- "N"; lab <- sapply(strsplit(lab, char), function (x) {if (length(x)>1) {gsub(" ","",paste(x[[1]], char, collapse = ""))} else {x}})
#cols <- c("purple","green","royalblue","darkorange","darkgreen","red","black","magenta","brown","forestgreen","darkblue","blue")
#text(x[s], y[s], labels = lab, cex=1.3, col = cols[as.numeric(annot[s,]$Donor)])

# TODO Kernel-PCA
#library(kernlab)
#kpc <- kpca(~., data=as.data.frame(x.train), kernel="rbfdot", features = 2)
#Print the principal component vectors
#pcv(kpc)
#Plot the data projection on the components
#par(mfrow=c(2,2))
#plot(rotated(kpc), col = as.integer(as.factor(ann[rownames(x.train),"time"])), xlab="1st Principal Component", ylab="2nd Principal Component")
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



#' Description: Barplot of top findings from pairwise.comparisons function
#'              
#' Arguments:
#'   @param top.findings Ouput from pairwise.comparisons function
#'   @param topN number of top findings to plot
#'   @param annot annotation matrix    
#'   @param smat species abundancy matrix
#' Returns:
#'   @return List with top findings from pairwise comparisons and their q-values
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

top.barplots <- function (top.findings, topN = 5, annot, smat) {

  # Investigate top findings
  specs <- unique(unlist(sapply(top.findings, rownames)))[1:topN] # phylotypes

  if (length(specs)>0) {
    specs.shortnames <- sapply(specs, function(x) {ss <- strsplit(x, " "); if (length(ss[[1]]) > 1) { paste(paste(substr(ss[[1]][[1]], 1, 1), ".", sep = ""), paste(ss[[1]][-1], collapse = " "), sep = " ")} else {ss[[1]][[1]]}})
    names(specs.shortnames) <- specs
    s <- rownames(smat) # samples
    df <- cbind(smat[s, specs], annot[s,])
    means <- reshape::melt(aggregate(df[specs], by=list(varname = df[[varname]]), FUN=mean))
    stds <- reshape::melt(aggregate(df[specs], by=list(varname = df[[varname]]), FUN=sd))
    Nsqrt <- sqrt(as.numeric(reshape::melt(aggregate(df[specs], by = list(varname = df$varname), FUN = length))$value))
    dfm <- cbind(means, mean = means$value, sd = stds$value, sd.of.mean = as.numeric(stds$value) / Nsqrt)
    dfm$shortnames <- as.factor(specs.shortnames[dfm$variable])

    # Create the barplot component
    p <- ggplot2::ggplot(dfm, aes(x = shortnames, y = value, fill = varname)) 
    dodge <- position_dodge(width = 0.9)
    p <- p + geom_bar(position="dodge") 
    p <- p + geom_errorbar(aes(x = shortnames, ymax = mean + 1.96*sd.of.mean, ymin = mean - 1.96*sd.of.mean), position = dodge, width=0.25)
    p <- p + scale_fill_grey() + theme_bw() + ylab("Signal") + xlab("") 
    p <- p + opts(axis.text.x=theme_text(angle=-20)) 

    } else {
      warning("No features available.")
      p <- NULL  
    }

}




#' Description: Visualize top findings from pairwise.comparisons 
#'          
#' FIXME: merge with top.barplots?
#'     
#' Arguments:
#'   @param x data matrix
#'   @param y annotaction factor
#'   @param oligo.map mapping between features
#'   @param color.level feature level to color
#'   @param bar.level feature level for bars
#'   @param top.findings output from pairwise.comparisons function
#'   @param top.findings.qvals output from pairwise.comparisons function
#'   @param qth q-value threshold
#'   @param qth.star q-value threshold for stars
#'   @param mode barplot / heatmap
#' Returns:
#'   @return List with top findings from pairwise comparisons and their q-values
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PlotTopComparisons <- function (x, y, oligo.map, color.level, bar.level, top.findings, top.findings.qvals, qth = 0.05, qth.star = 0.01, mode = "barplot") {

  # x <- smat; color.level = "level 1"; bar.level = "level 2"; qth <- 0.05; mode = "heatmap"

  flevels <- as.character(unique(y))

  map <- oligo.map[[color.level]]
  names(map) <- as.character(oligo.map[[bar.level]])

  # Pick the most significant findings only
  nams <- unique(unlist(sapply(top.findings, function(tab){unlist(rownames(subset(tab, qvalue < qth)))})))

  x <- x[, nams]  

  fc <- aggregate(x, by = list(varname = y), mean) # Mean for each group

  # Compare each group to other group
  comparisons <- list()

  for (i in 1:(length(flevels)-1)) {
    for (j in (i+1):length(flevels)) {
      di <- flevels[[i]]
      dj <- flevels[[j]]
      nam <- paste(di, dj, sep = "-")
      comparisons[[nam]] <- as.numeric(subset(fc, varname == di)[-1]) - as.numeric(subset(fc, varname == dj)[-1])
    }
  }

  nams <- names(comparisons)
  df <- data.frame(comparisons)
  colnames(df) <- nams
  rownames(df) <- colnames(x)
  df$bar.level <- colnames(x)
  df$color.level <- factor(map[as.character(df$bar.level)]) # Check L2 for the phylotype species

  if (mode == "barplot") {

    specs.shortnames <- sapply(as.character(df$bar.level), function(x) {ss <- strsplit(x, " "); if (length(ss[[1]]) > 1) { paste(paste(substr(ss[[1]][[1]], 1, 1), ".", sep = ""), paste(ss[[1]][-1], collapse = " "), sep = " ")} else {ss[[1]][[1]]}})
    names(specs.shortnames) <- as.character(df$bar.level)
    df$shortnames <- factor(specs.shortnames[df$bar.level])
    dfm <- reshape::melt(df)

    p <- ggplot2::ggplot(dfm) 
    dodge <- position_dodge(width = 0.9) # define the width of the dodge
    p <- p + geom_bar(position = "dodge", width = 0.4, aes(x = shortnames, y = value, fill = color.level))
    p <- p + facet_grid(.~variable, scales = "free", space = "free")
    p <- p + coord_flip()
    p <- p + theme_bw() + ylab("") + xlab("Fold Change")
    p <- p + opts(axis.text.x = theme_text(size = 7))
    # p <- p + opts(legend.position = c(0.15,.15))

  } else if (mode == "heatmap") {

    ncomp <- ncol(df) - 2 # number of comparisons

    dfm <- reshape::melt(df[, 1:(ncomp + 1)], variable = "comparison") # reshape::melt matrix
    keep <- !is.na(dfm$value)
    dfm <- dfm[keep, ]
    dfm$comparison <- droplevels(dfm$comparison)
    #flevels <- levels(dfm$comparison)

    qmat <- NULL; nams <- c()
    if (length(flevels)>1) {
      for (i in 1:(length(flevels)-1)) {
        for (j in (i+1):(length(flevels))) {
  
          ai <- flevels[[i]]
          aj <- flevels[[j]]
          nam <- paste(ai, aj, sep = "-")

          qmat <- cbind(qmat, unlist(top.findings.qvals[[ai]][[aj]])[df[[bar.level]]])
          nams <- c(nams, nam)
        }
      }
    } else {

          nam <- flevels[[1]]
          qmat <- top.findings.qvals[[1]][[1]][df$bar.level]
	  qmat <- matrix(qmat, length(qmat))
	  rownames(qmat) <- df$bar.level
          nams <- c(nams, nam)

    }

    if (nrow(qmat) > 0) {

      qmat <- data.frame(qmat)
      colnames(qmat) <- nams
 
      qmat[["phylotype"]] <- rownames(qmat)
      dfm.qval <- reshape::melt(qmat)
      names(dfm.qval) <- "value"

      cex.xlab = 7; cex.ylab = 12

      p <- ggplot2::ggplot(dfm, aes(x = comparison, y = bar.level, fill = value))
      p <- p + geom_tile() 
      limit <- ceiling(max(abs(na.omit(dfm$value))))
      p <- p + scale_fill_gradientn("Fold Change", breaks=seq(from=limit, to=-limit, by=-0.2), colours = c("darkblue", "blue", "white", "red", "darkred"), limits=c(-limit,limit)) 
      p <- p + theme_bw() 
      p <- p + opts(axis.text.x=theme_text(angle = 0, size = cex.xlab)) 
      p <- p + opts(axis.text.y=theme_text(size = cex.ylab))

      # Merkkaa merkitsevat tahdella
      p <- p + geom_text(data = subset(dfm.qval, value < qth.star), aes(x = variable, y = value, label = "+"), col = "white", size = 3)
      p <- p + xlab("") + ylab("")
    }
  }

  p

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




#' Description: Calculate and plot hierarchical clustering for profiling 
#' script output
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param data.dir data directory
#'   @param phylogeny phylogeny
#'
#' Returns:
#'   @return hclust object
#'
#' @export
#' @examples # dat <- read.profiling(params$wdir, "species", "rpa"); hc <- add.hclust.plots(dat)
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

add.hclust.plots <- function (dat, data.dir, phylogeny) {

  # Read heatmap plotting parameters		 
  hc.params <- ReadHclustParameters(dat, data.dir)

  # Plot and save into output file
  message(paste("Storing hclust image in", hc.params$file))  
  plotdev <- png(filename = hc.params[["file"]], width = max(trunc(hc.params[["ppcm"]]*21), trunc(hc.params[["ppcm"]]*21*ncol(dat)/70)), height=trunc(hc.params[["ppcm"]]*29.7)) 
  try(PlotPhylochipHeatmap(data = dat,
                phylogeny = phylogeny,
                metric = hc.params[["clmet"]],
                tax.level = hc.params[["lev"]],
                include.tree = ifelse(hc.params[["tree.display"]] == 'yes', TRUE, FALSE),
                palette = hc.params[["pal"]],
                fontsize = hc.params[["fontsize"]],
                figureratio = hc.params[["figureratio"]], 
		hclust.method = hc.params[["hclust.method"]])) 
  dev.off()


  # Plot CLUSTER TREE TO A GRAPHICS WINDOW
  if(ncol(dat) > 2){
  
    hc <- calculate.hclust(log10(dat + 1), hclust.method = hc.params[["hclust.method"]], metric = list.clustering.metrics()[[hc.params[["clmet"]]]])

    x11() 
    plot(hc$log10, hang = -1, main = "Hierarchical clustering, oligolevel, log10")

    x11() 
    plot(hc$raw, hang = -1, main = "Hierarchical clustering, oligolevel, raw")

  } else {

    warning("Clustering skipped because there are too few samples.\n")

  }

  list(ppcm = ppcm, hclust.method = hclust.method, pal = pal, lev = lev, clmet = clmet, tree.display = tree.display, figureratio = figureratio, fontsize = fontsize)

}

