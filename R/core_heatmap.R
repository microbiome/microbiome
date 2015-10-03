#' Core heatmap
#'
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param detection.thresholds A vector or  a scalar indicating the number of intervals
#'        in (0, log10(max(data))). The detection thresholds are calculated for relative abundancies.
#' @param palette palette for the plot.type = 'heatmap'
#'  
#' @return Used for its side effects
#' @importFrom RColorBrewer brewer.pal
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
core_heatmap <- function(x, detection.thresholds = 20, palette = "bw") {

    # Calculate relative abundances
    x <- transform_sample_counts(x, function(x) x/sum(x))

    # Get OTU matrix
    data <- t(otu_table(x)@.Data)

    DetectionThreshold <- Taxa <- Prevalence <- NULL
    
    if (length(detection.thresholds) == 1) {
      detection.thresholds <- 10^seq(log10(1e-3), log10(max(data)), length = detection.thresholds)
   }
    
    # Prevalences with varying detection thresholds
    prev <- lapply(detection.thresholds, function (th) {prevalence(data, detection.threshold = th)})
    prev <- 100*do.call("cbind", prev)
    colnames(prev) <- as.character(detection.thresholds)

    df <- melt(prev)
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    o <- names(sort(rowSums(prev)))
    df$Taxa <- factor(df$Taxa, levels = o)
    df$DetectionThreshold <- 100*df$DetectionThreshold

    theme_set(theme_bw(10))
    p <- ggplot(df, aes(x = DetectionThreshold, y = Taxa, fill = Prevalence))
    p <- p + scale_x_log10()
    p <- p + geom_tile()
    p <- p + xlab("Detection Threshold (Relative Abundance %)")

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
