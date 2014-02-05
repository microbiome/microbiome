#' add.heatmap
#' Description: Add oligprofile heatmap into output directory
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param output.dir output data directory
#'   @param output.file output file name
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
#' @examples # data(peerj32); hc <- add.heatmap(peerj32$microbes[, 1:4])
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

add.heatmap <- function (dat, output.dir, output.file = NULL, phylogeny.info, ppcm = 150, 
	         hclust.method = "complete", palette = "white/black", level = "L1", metric = "pearson", 
  		 figureratio = 10, fontsize = 40, tree.display = TRUE) {

  # dat <- finaldata[["oligo"]]; output.dir = params$wdir;  output.file = NULL; phylogeny.info = phylogeny.info; ppcm = 150; hclust.method = "complete"; palette = "white/blue"; level = "L2"; metric = "pearson"; figureratio = 12; fontsize = 12; tree.display = TRUE
  #output.dir = "~/tmp/";  output.file = NULL; phylogeny.info = phylogeny.info; ppcm = 150; hclust.method = "complete"; palette = "white/blue"; level = "L2"; metric = "pearson"; figureratio = 12; fontsize = 12; tree.display = TRUE

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
#'   @param hclust.method hierarchical clustering method. See help(hclust) for details. To prevent ordering of the rows, use hclust.method = NULL.
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

  if (is.null(hclust.method) || hclust.method == "none") {hclust.method <- NULL}

  params <- c(metric = metric, 
  	      level = level, 
	      tree.display = tree.display, 
	      palette = palette, 
  	      fontsize = fontsize, 
	      figureratio = figureratio,  
	      hclust.method = hclust.method)
			 
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
   if (metric == "euclidean" && !is.null(hclust.method)) {
     hc <- hclust(dist(t(data)), method = hclust.method)
     ord <- hc$order
   } else if (metric == "pearson" && !is.null(hclust.method)) {
     hc <- hclust(as.dist(1 - cor(data, use = "pairwise.complete.obs")), method = hclust.method)
     ord <- hc$order
   } else if (metric == "spearman" && !is.null(hclust.method)) {
     hc <- hclust(as.dist(1 - cor(data, use = "pairwise.complete.obs", method = "spearman")), method = hclust.method)
     ord <- hc$order
   } else if (is.null(hclust.method)) {
     ord <- 1:ncol(data) # do not change the order
   }

   # Order the data
   data <- data[, ord] 

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

   if (tree.display && !is.null(hclust.method)) {
      layout(matrix(c(3,0,1,2), ncol = 2, byrow = TRUE),widths=lcm(c(profilewidth,annwidth)),heights=lcm(figureheight*heights))
   } else {    
      layout(matrix(c(3,0,1,2),ncol=2,byrow=TRUE),widths=lcm(c(profilewidth,annwidth)),heights=lcm(figureheight*heights))
   }

   par(mar = c(1,1,0,0), oma = c(0,0,0,0))
   
   img <- as.matrix(rev(as.data.frame(t(data[as.character(phylogeny.info$oligoID),]))))
   image(z = img, col = palette, axes = FALSE, frame.plot = TRUE, zlim = limits)
   plot.new()
   par(mar = c(1, 0, 0, 1), usr = c(0, 1, 0, 1), xaxs = 'i', yaxs = 'i')

   rect(xleft = rep(0,length(levs)),
        ybottom = c(0,cumsum(rev(levs))[1:length(levs)-1])/sum(levs),
        xright = rep(1,length(levs)),ytop=cumsum(rev(levs))/sum(levs), border = 'grey')

   text(x = c(0.03), y = (cumsum(rev(levs))-rev(levs/2))/sum(levs), labels = (names(rev(levs))), pos = 4)

   if (tree.display && !is.null(hclust.method)) {

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


#' Description: List color scales
#'
#' Arguments:
#'
#' Returns:
#'   @return list of color scales
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal

list.color.scales <- function () {
  ## Different colour scales
  list('white/blue'=colorRampPalette(c("white","darkblue"),interpolate='linear')(100),
       'white/black'=colorRampPalette(c("white","black"),interpolate='linear')(100),
       'black/yellow/white'=colorRampPalette(c("black","yellow","white"),bias=0.5,interpolate='linear')(100))

}


