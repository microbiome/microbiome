#' @title OTU core heatmap
#' @description Core heatmap.
#' @param data OTU matrix
#' @param detection.thresholds A vector or a scalar indicating the number of intervals in (0, log10(max(data))). The detection thresholds are calculated for relative abundancies.
#' @param palette palette for the plot.type = 'heatmap'
#' @param min.prevalence If minimum prevalence is set, then filter out those rows (taxa) and columns (detection thresholds) that never exceed this prevalence threshold. This helps to zoom in on the actual core region of the heatmap.
#' @param taxa.order Ordering of the taxa.
#' @return Used for its side effects
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_heatmap <- function(data, detection.thresholds = 20, palette = "gray", min.prevalence = NULL, taxa.order = NULL, ncolor = 5) {

    DetectionThreshold <- Taxa <- Prevalence <- NULL
    
    # Prevalences with varying detection thresholds
    prev <- lapply(detection.thresholds, function (th) {prevalence(data, detection.threshold = th)})
    prev <- do.call("cbind", prev)
    colnames(prev) <- as.character(detection.thresholds)

    # Exclude rows and cols that never exceed the given prevalence
    if (!is.null(min.prevalence)) {
      prev <- prev[rowMeans(prev > min.prevalence) > 0, colMeans(prev > min.prevalence) > 0]
    }
    
    df <- as.data.frame(prev)
    df$ID <- rownames(prev)
    df <- gather(df, "ID")
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    df$DetectionThreshold <- as.numeric(as.character(df$DetectionThreshold))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))

    if (is.null(taxa.order)) {
      o <- names(sort(rowSums(prev)))
    } else {
      o <- taxa.order
    }
    df$Taxa <- factor(df$Taxa, levels = o)

    theme_set(theme_bw(10))
    p <- ggplot(df, aes(x = DetectionThreshold, y = Taxa, fill = Prevalence))
    p <- p + geom_tile()
    p <- p + xlab("Detection Threshold")    
    p <- p + scale_x_log10()


    colours = palette(seq(0, 1, length = ncolor))
    #if (palette == "bw") {
    #  #colours <- c("black", "darkgray", "gray", "lightgray", "white")
    #} else if (palette == "spectral") {
    #  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    #  colours <- myPalette(5)
    #}
    
    p <- p + scale_fill_gradientn("Prevalence", 
        breaks = seq(from = 0, to = 100, by = 10),
	colours = colours, limits = c(0, 100))
    
    return(list(plot = p, data = df))
    
} 

