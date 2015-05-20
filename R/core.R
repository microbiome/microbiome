#' core_matrix
#'
#' create core matrix 
#'
#' @param x \code{\link{phyloseq}} object
#' @param verbose verbose
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range
#'
#' @return Estimated core microbiota
#'
#' @examples 
#'   #pseq <- download_microbiome("peerj32")$physeq
#'   #core <- core_matrix(pseq)
#'
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_matrix <- function(x,  
          prevalence.intervals = seq(5, 100, 5), 
          detection.thresholds = NULL, verbose = FALSE) {
    
    # Convert into OTU matrix
    data <- otu_table(x)@.Data

    # Use log10 
    data <- log10(data)

    # Convert prevalences from percentages to numerics
    p.seq <- 0.01 * prevalence.intervals * ncol(data)

    ## Intensity vector
    if (is.null(detection.thresholds)) {
      i.seq <- seq(min(data), max(data), length = 10)
    } else {   
      # Use log10
      i.seq <- log10(detection.thresholds)
    }
    
    coreMat <- matrix(NA, nrow = length(i.seq), ncol = length(p.seq), 
                      	  dimnames = list(i.seq, p.seq))
    
    n <- length(i.seq) * length(p.seq)
    cnt <- 0
    for (i in i.seq) {
      for (p in p.seq) {
        if (verbose) {
          cnt <- cnt + 1
        }
          coreMat[as.character(i), as.character(p)] <- core.sum(data, i, p)
        }
    }
    
    # Convert Prevalences to percentages
    colnames(coreMat) <- 100 * as.numeric(colnames(coreMat))/ncol(data)
    rownames(coreMat) <- as.character(10^as.numeric(rownames(coreMat)))
    
    coreMat

}

core.sum <- function(data, intTr, prevalenceTr) {

    d.bin <- data > intTr
    prevalences <- rowSums(d.bin)
    nOTUs <- sum(prevalences >= prevalenceTr)
    return(nOTUs)

}



#' plot_core
#'
#' Core visualization 2D
#'
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param title title
#' @param plot plot the figure 
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range
#' @param plot.type Plot type ('lineplot' or 'heatmap')
#' @param palette palette for the plot.type = 'heatmap'
#'  
#' @return Used for its side effects
#'
#' @examples 
#' #pseq <- download_microbiome("atlas1006")
#' #p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)))
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_core <- function(x, title = "Core", plot = TRUE, 
		   prevalence.intervals = NULL, 
		   detection.thresholds = NULL, 
		   plot.type = "lineplot", palette = "bw") {
    
  if (plot.type == "lineplot") {
    p <- core_lineplot(x,  
      	 	       prevalence.intervals = prevalence.intervals, 
		       detection.thresholds = detection.thresholds) 

  } else if (plot.type == "heatmap") {
    res <- core_heatmap(x, detection.thresholds = detection.thresholds, palette = palette)
    p <- res$p
    # Data is available but not returned in current implementation
    prevalences <- res$prevalences
  }

  p <- p + ggtitle(title)
    
  if (plot) {
    print(p)
  }

  p
}


core_lineplot <- function (x, title = "Common core",  
                   xlabel = "Abundance", 
                   ylabel = "Core size (number of taxa)", 
		   prevalence.intervals = NULL, 
		   detection.thresholds = NULL) {

    if (class(x) == "phyloseq") {
      x <- core_matrix(x, prevalence.intervals, detection.thresholds)
    }

    Abundance <- Prevalence <- Count <- NULL
    
    df <- melt(x)
    names(df) <- c("Abundance", "Prevalence", "Count")

    theme_set(theme_bw(20))
    p <- ggplot(df, aes(x = Abundance, y = Count, color = Prevalence, 
                        group = Prevalence))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + scale_x_log10()
    p <- p + xlab(xlabel)
    p <- p + ylab(ylabel)
    p <- p + ggtitle(title)
        
    return(p)
}


#' core_heatmap
#'
#' Core heatmap
#'
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param detection.thresholds a vector of intensities around the data range
#' @param palette palette for the plot.type = 'heatmap'
#'  
#' @return Used for its side effects
#' @importFrom RColorBrewer brewer.pal
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_heatmap <- function(x, detection.thresholds = NULL, palette = "bw") {

    # Convert into OTU matrix
    data <- t(otu_table(x)@.Data)

    DetectionThreshold <- Taxa <- Prevalence <- NULL
    
    if (is.null(detection.thresholds)) {
      detection.thresholds <- seq(min(data), max(data), length = 10)
    }
    
    # Prevalences with varying detection thresholds
    prev <- lapply(detection.thresholds, function (th) {prevalence(data, detection.threshold = th)})
    prev <- 100*do.call("cbind", prev)
    colnames(prev) <- as.character(detection.thresholds)

    df <- melt(prev)
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    o <- names(sort(rowSums(prev)))
    df$Taxa <- factor(df$Taxa, levels = o)

    theme_set(theme_bw(10))
    p <- ggplot(df, aes(x = DetectionThreshold, y = Taxa, fill = Prevalence))
    p <- p + scale_x_log10()
    p <- p + geom_tile()
    
    if (palette == "bw") {
        colours <- c("black", "darkgray", "gray", "lightgray", "white")
    } else if (palette == "spectral") {
        myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
        colours <- myPalette(5)
    }
    
    p <- p + scale_fill_gradientn("Prevalence", 
        breaks = seq(from = 0, to = 100, 
        by = 10), colours = colours, limits = c(0, 100))
    
    return(list(plot = p, data = prev))
    
} 
