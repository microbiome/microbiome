#' Description: Calculate and plot hierarchical clustering for profiling script output
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param data.dir data directory
#'   @param phylogeny phylogeny
#' Returns:
#'   @return hclust object
#'
#' @export
#' @examples # dat <- read.profiling(params$wdir, "species", "rpa"); hc <- add.hclust.plots(dat)
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
    plot(hc$log10, hang=-1, main="Hierarchical clustering, oligolevel, log10")

    x11() 
    plot(hc$raw, hang=-1, main="Hierarchical clustering, oligolevel, raw")

  } else {

    warning("Clustering skipped because there are too few samples.\n")

  }

  list(ppcm = ppcm, hclust.method = hclust.method, pal = pal, lev = lev, clmet = clmet, tree.display = tree.display, figureratio = figureratio, fontsize = fontsize)

}




#' Description: Read parameters for PROFILE PLOT FUNCTION
#'
#' Arguments:
#'   @param dat data matrix
#'   @param data.dir data directory
#'   @param ppcm plotting parameter
#'   @param hclust.method hierarchical clustering method
#'   @param pal color palette
#'   @param lev phylogeny level
#'   @param clmet clustering metrics
#'   @param tree.display tree.display
#'   @param figureratio figureratio
#'   @param fontsize fontsize
#'
#' Returns:
#'   @return plots the profile
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

ReadHclustParameters <- function (dat, data.dir, ppcm = 150, hclust.method = "complete", pal = "white/blue", lev = "level 2", clmet = "Pearsons correlation coefficient", tree.display = "yes", figureratio = 12, fontsize = 12) {

  cmetrics <- list.clustering.metrics()
  cscales  <- list.color.scales()
  include.tree <- TRUE

  ## Plotting parameters
  defParams <- tk_select.list(c("Yes", "No"), preselect="Yes", multiple = FALSE,
                              title = "Use default plotting parameters?")

   if(defParams=="No"){

     pal <- tk_select.list(names(cscales),preselect=names(cscales)[1],
                         multiple=FALSE,title="Select colour scale")

     lev <- tk_select.list(names(dat), preselect = names(dat)[[1]],
                         multiple = FALSE, title = "Select taxonomic level for the graph")
   
     clmet <- tk_select.list(names(cmetrics),preselect=names(cmetrics)[1],
                           multiple=FALSE,title="Select clustering metrics")

     hclust.method <- tk_select.list(c("complete", "ward"),preselect = "complete",
                           multiple=FALSE,title="Select clustering method")

     tree.display <- tk_select.list(c("yes","no"),preselect=c("yes"),
                                  multiple=FALSE,title="Display tree in graph?")

     figureratio <- NA
     if (tree.display == "yes") {
       figureratio <- tk_select.list(c(15,12,10,8,6,4),
					preselect = 12,
                                   	multiple = FALSE,
				   	title = "Percentage of tree height to total figure")
     }

     fontsize <- tk_select.list(c(8:12, seq(14, 40, 2)),preselect = '9',multiple = FALSE, title = "Font size")

   } 

  # Graph to be saved into 
  clusterGraphFile <- paste(data.dir,"/", gsub(" ", "", lev), "-oligoprofileClustering.png",sep="")

  ## make the figure width as a function of the number of the samples
  if(ncol(dat) < 3 ) { tree.display <- 'no' }

  list(file = clusterGraphFile, ppcm = ppcm, clmet = clmet, lev = lev, include.tree = include.tree, palette = pal, fontsize = fontsize, figureratio = figureratio, hclust.method = hclust.method, tree.display = tree.display)

}

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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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

#' Description: transparent color
#'
#' Arguments:
#'   @param color color
#'   @param tr transparency level
#'
#' Returns:
#'   @return TBA
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

trans <- function(color, tr = 150){
  trCol <- rgb(col2rgb(color)[1],
               col2rgb(color)[2],
               col2rgb(color)[3],tr, maxColorValue=255)
}


#' Description: Plot PCA with labels and explained variances
#'
#' Arguments:
#'   @param dPCA TBA
#'   @param main TBA
#'   @param xPC TBA
#'   @param yPC TBA
#'   @param showLabels TBA
#'   @param type TBA
#'   @param ... further parameters to be passed
#'
#' Returns:
#'   @return TBA
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

plotPCA <- function(dPCA, main="PCA plot", xPC = 1, yPC = 2, showLabels = T, type = "n", ...){
   x11(width=11, height=11)
   dPCA.propVars <- round(dPCA$sdev^2/sum(dPCA$sdev^2)*100,digits=2)
   plot(dPCA$x[,xPC],dPCA$x[,yPC],type=type, main=main,
        xlab=paste("PC ",xPC," (",dPCA.propVars[xPC],"% of variance)",sep=""),
        ylab=paste("PC ",yPC," (",dPCA.propVars[yPC],"% of variance)",sep=""),...)
   if(showLabels)
     text(dPCA$x[,xPC],dPCA$x[,yPC]-0.02*max(dPCA$x[,yPC]),
          labels=rownames(dPCA$x), cex=0.7,...)
}

#' Description: 
#'
#' Arguments:
#'   @param dPCA TBA
#'   @param groups TBA
#'   @param main TBA
#'   @param xPC TBA
#'   @param yPC TBA
#'   @param type TBA
#'   @param legLoc TBA
#'   @param w TBA
#'   @param h TBA
#'   @param ... further parameters to be passed
#'
#' Returns:
#'   @return plot
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

plotGroupPCA <- function(dPCA,groups,main="PCA plot",xPC=1, yPC=2,type="n",
                         legLoc="topright",w=11,h=11,...){
   x11(width=w, height=h)
   dPCA.propVars <- round(dPCA$sdev^2/sum(dPCA$sdev^2)*100,digits=2)
   plot(dPCA$x[,xPC],dPCA$x[,yPC],type=type, main=main,
        xlab=paste("PC ",xPC," (",dPCA.propVars[xPC],"% of variance)",sep=""),
        ylab=paste("PC ",yPC," (",dPCA.propVars[yPC],"% of variance)",sep=""),...)
   ##text(dPCA$x[,xPC],dPCA$x[,yPC]-0.02*max(dPCA$x[,yPC]),
   ##    labels=rownames(dPCA$x), cex=0.7,...)
   uGroups <- unique(groups)
   nG <- length(uGroups)
   gCols <- rainbow(nG)
   for(g in 1:nG){
     gFilter <- groups==uGroups[g]
     points(dPCA$x[gFilter,xPC],dPCA$x[gFilter,yPC],col=gCols[g], pch=19) 
   }
   legend(legLoc,legend=uGroups, fill=gCols)
 }
 
 

#' Description: Plot oligoprofile figures without data retrieval from the db
#'
#' Arguments:
#'   @param d data
#'   @param tax.level Taxonomic level
#'   @param metric metrics
#'   @param include.tree TBA
#'   @param fontsize font size
#'   @param figureratio figure ratio
#'   @param oligomap oligomap
#'   @param width width
#'   @param height height
#'   @param oligomapJN TBA
#'   @param plot.profile TBA
#'
#' Returns:
#'   @return plot
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

plotOligoP <- function(d, tax.level='level1',
                        metric='correlation', include.tree=T,
                        fontsize=11, figureratio=15, oligomap=oligomapJN,
                        width=12, height=20, oligomapJN = NULL, 
			plot.profile = NULL){
 
   par(mar=c(0.1,0.1,0.1,0.1))
   x11(width=width, height=height)
   plot.profile(data=log(d), metric=metric,
                oligomap=oligomap, tax.level=tax.level,
                fontsize=fontsize, include.tree=include.tree,
                figureratio=figureratio)
}


#' Description: Plot Mean Profile 
#'
#' Arguments:
#'   @param data TBA
#'
#' Returns:
#'   @return barplot
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

plotMeanProfile <- function(data){
  require(gplots)
  data.mean <- apply(data,1, mean)
  data.sd <- apply(data, 1, sd)
  barplot2(data.mean, horiz=T,
           plot.ci=T, ci.u=data.mean+data.sd, ci.l=data.mean-data.sd)
}


#' Description: plot scatter
#'
#' Arguments:
#'   @param sampleA TBA
#'   @param sampleB TBA
#'   @param d TBA
#'   @param cex TBA
#'
#' Returns:
#'   @return scatterplot 
#'
#' @export 
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

plotScatter <- function(sampleA, sampleB, d, cex=1){
  x11()
  plot(d[,sampleA],
       d[,sampleB], pch=".",
       xlab=sampleA, ylab=sampleB,
       main=paste("Between sample similarity, Pearson: ",
         round(cor(d[,sampleA],d[,sampleB]),2)))
  points(d[,sampleA],
         d[,sampleB], cex=cex, pch=".")
}



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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
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

