# Copyright (C) 2011-2014 Jarkko Salojarvi and Leo Lahti 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: core.sum
#'
#' Arguments:
#'   @param data data matrix; phylotypes vs. samples
#'   @param intTr intTr
#'   @param prevalenceTr prevalenceTr
#'
#' Returns:
#'   @return TBA
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

core.sum <- function(data, intTr, prevalenceTr){
  d.bin <- data>intTr
  prevalences <- rowSums(d.bin)
  # jos haluat tietaa lajit, ala summaa!
  nOTUs <- sum(prevalences>=prevalenceTr) 
  return(nOTUs)
}


#' Description: core.which
#'
#' Arguments:
#'   @param data data matrix; phylotypes vs. samples
#'   @param intTr intTr
#'   @param prevalenceTr prevalenceTr
#'
#' Returns:
#'   @return TBA
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

core.which <- function(data, intTr, prevalenceTr){
  d.bin <- data>=intTr
  prevalences <- rowSums(d.bin)
  nOTUs <- as.numeric(prevalences>=prevalenceTr)
  return(nOTUs)
}



#' createCore
#'
#' Description: create coreMatrix 
#'
#' Arguments:
#'   @param data data matrix; phylotypes vs. samples
#'   @param verbose verbose
#'   @param prevalence.intervals a vector of prevalence percentages in [0,100]
#'   @param intensity.intervals a vector of intensities around the data range
#'
#' Returns:
#'   @return TBA
#'
#' @examples data(peerj32); 
#'          core <- createCore(t(peerj32$microbes))
#'
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

createCore <- function(data, verbose = FALSE, 
                   prevalence.intervals = seq(20, 100, 20), 
               intensity.intervals = NULL) {

   ## Prevalence vector
   if (is.null(prevalence.intervals)) {
     prevalence.intervals <- seq(0, 100, 10)
   }
   # Convert prevalences from percentages to numerics
   p.seq <- 0.01*prevalence.intervals*ncol(data)

   ## Intensity vector
   if (is.null(intensity.intervals)) {
     i.seq <- seq(min(data), max(data), length = 10) 
   } else {
     i.seq <- intensity.intervals
   }

   coreMat <- matrix(NA, nrow = length(i.seq), ncol = length(p.seq), 
                      dimnames = list(i.seq, p.seq))

   n <- length(i.seq)*length(p.seq)
   cnt <- 0
   for(i in i.seq){
     for(p in p.seq){
       if (verbose) {cnt <- cnt + 1; message(cnt/n)}
       coreMat[as.character(i),as.character(p)] <- core.sum(data, i, p)
     }
   }

   # Convert Prevalences to percentages
   colnames(coreMat) <- 100*as.numeric(colnames(coreMat))/ncol(data)

   return(coreMat)
}



#' Core3D
#'
#' Description: Core visualization 3D
#'
#' Arguments:
#'  @param coreMat core matrix
#'  @param title title
#'  @param xlab X axis label
#'  @param cex.axis axis text size
#'
#' Returns:
#'  @return Used for its side effects
#'
#' @examples data(peerj32); 
#'          c3d <- Core3D(createCore(t(peerj32$microbes)))
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Core3D <- function(coreMat, title = "Core microbiota", 
                         xlab = "Minimum Intensity", cex.axis = 0.7){

  MinimumPrevalence <- as.numeric(colnames(coreMat))
  MinimumLogIntensity <- as.numeric(rownames(coreMat))
  tmp <- persp(MinimumLogIntensity, MinimumPrevalence, coreMat, 
                 theta = 60, phi = 5, 
           main=title, 
           col="light blue", 
           axes=T, ticktype="detailed", nticks=9, 
           shade = 0.58,         
           cex.axis = cex.axis, 
           ylab = "Minimum Prevalence", xlab = xlab, zlab = "Core Size")

  return(NULL)

}


#' Core2D
#'
#' Description: Core visualization 2D
#'
#' Arguments:
#'   @param coreMat core matrix
#'   @param title title
#'   @param plot plot the figure 
#'   @param xlabel X axis label
#'   @param ylabel Y axis label
#'  
#'
#' Returns:
#'   @return Used for its side effects
#'
#' @examples data(peerj32); 
#'          c2d <- Core2D(createCore(t(peerj32$microbes)))
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Core2D <- function(coreMat, title = "Common core", plot = TRUE, 
                      xlabel = "Abundance", ylabel = "Core size (number of taxa)"){

  Abundance <- Prevalence <- Count <- NULL       

  df <- melt(coreMat)
  names(df) <- c("Abundance", "Prevalence", "Count")
  theme_set(theme_bw(20))
  p <- ggplot(df, aes(x = Abundance, y = Count, color = Prevalence, 
                 group = Prevalence))
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + xlab(xlabel)
  p <- p + ylab(ylabel)
  p <- p + ggtitle("Core microbiota")

  if (plot) { print(p) }

  return(p)

}


#' core_heatmap
#'
#' Description: Heatmap of core microbiota
#'
#' Arguments:
#'   @param data data matrix: phylotypes vs. samples
#'   @param detection.thresholds Vector of detection thresholds
#'   @param plot plot the figure
#'   @param palette "bw" (grayscale) or "spectral" (colourscale)
#'
#' Returns:
#'   @return List with the following elements: 
#'            plot: ggplot figure
#'         data: prevalence data with the varying thresholds
#'
#' @examples data(peerj32); 
#'          core <- core_heatmap(t(peerj32$microbes))
#'
#' @export 
#' @importFrom reshape melt
#' @import ggplot2
#' @import RColorBrewer
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

core_heatmap <- function (data, detection.thresholds = NULL, 
                              plot = TRUE, 
                  palette = "bw") {

  DetectionThreshold <- Taxa <- Prevalence <- NULL

  if (is.null(detection.thresholds)) {
    detection.thresholds <- seq(min(data), max(data), length = 10)
  }

  # Prevalences with varying detection thresholds
  taxa <- rownames(data)
  prevalences <- matrix(NA, nrow = length(taxa), 
            ncol = length(detection.thresholds))
  rownames(prevalences) <- taxa
  colnames(prevalences) <- as.character(detection.thresholds)
  for (det.th in detection.thresholds) {
    prevalence <- 100*sort(rowMeans(data > det.th))
    prevalences[taxa, as.character(det.th)] <- prevalence[taxa]
  }

  df <- melt(prevalences)
  names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
  o <- names(sort(rowSums(prevalences)))
  df$Taxa <- factor(df$Taxa, levels = o)
  theme_set(theme_bw(10))
  p <- ggplot(df, aes(x = DetectionThreshold, y = Taxa, fill = Prevalence))
  p <- p + geom_tile()

  if (palette == "bw") {
    colours <- c("black", "darkgray", "gray", "lightgray", "white")
  } else if (palette == "spectral") {
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    colours <- myPalette(5)
  }

  p <- p + scale_fill_gradientn("Prevalence", 
                    breaks = seq(from = 0, to = 100, by = 10), 
            colours = colours,
            limits = c(0, 100))

  p <- p + ggtitle("Core microbiota") 

  if (plot) {
    print(p)
  }

  return(list(plot = p, data = prevalences))

}
