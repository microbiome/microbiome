#' htree.plot
#' 
#' Plot hierarchical clustering for the input data in absolute
#' and log10 scale using euclidean and pearson correlation similarities. 
#' Intended for internal use in the run.profiling.script function. 
#'
#' @param dat oligoprofile data in original (non-log) domain
#' @param method hierarchical clustering method
#' @param metric clustering similarity measure
#'
#' @return Used for its side effects; returns the arguments
#'
#' @export
#' @examples 
#'   data(peerj32)
#'   tmp <- htree.plot(peerj32$microbes[,1:5])
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

htree.plot <- function(dat, method = "complete", metric = "pearson") {
    
    # Plot CLUSTER TREES TO A GRAPHICS WINDOW Euclidean & Correlation
    #  / Raw & Log10
    
    if (ncol(dat) > 2) {
        
        if (metric == "euclidean") {
            # Metric: Euclidean
            hc <- hclust(dist(t(dat)), method = method)
            # hc.raw.eu <- hclust(dist(t(dat)), method = method) hc.log10.eu <-
            # hclust(dist(t(log10(dat + 1))), method = method)
        } else if (metric %in% c("pearson", "spearman")) {
            
            # 'Metric': Correlation
            hc <- hclust(as.dist(1 - cor(dat, use = "pairwise.complete.obs", 
                         method = metric)), 
                         method = method)
            
        }
        
        plot(hc, hang = -1, main = paste("hclust/", metric, sep = ""), 
             xlab = "Samples")
        
    } else {
      warning("Three or more samples required for clustering - skipped.\n")
    }
    
    list(data = dat, method = method, metric = metric)
    
}

