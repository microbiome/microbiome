#' plot_core
#'
#' Core visualization 2D
#'
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param title title
#' @param plot plot the figure 
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range,
#'          or a scalar indicating the number of intervals in the data range.
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
		   prevalence.intervals = seq(5, 100, 5), 		   
		   detection.thresholds = 20,
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
