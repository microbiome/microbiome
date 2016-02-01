#' @title Visualize OTU core
#' @description Core visualization 2D
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param title title
#' @param plot plot the figure 
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range,
#'          or a scalar indicating the number of intervals in the data range.
#' @param plot.type Plot type ('lineplot' or 'heatmap')
#' @param palette palette for the plot.type = 'heatmap'
#' @param min.prevalence If minimum prevalence is set, then filter out those rows (taxa) and columns (detection thresholds) that never exceed this prevalence threshold. This helps to zoom in on the actual core region of the heatmap. Only affects the plot.type = 'heatmap'.
#' @return A list with three elements: the ggplot object and the data. The data has a different form for the lineplot and heatmap. Finally, the applied parameters are returned.
#' @examples 
#' #pseq <- download_microbiome("atlas1006")
#' #p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)))
#' @export 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_core <- function(x, title = "Core", plot = TRUE, 
		   prevalence.intervals = seq(5, 100, 5), 		   
		   detection.thresholds = 20,
		   plot.type = "lineplot", palette = "bw", min.prevalence = NULL) {

  if (length(detection.thresholds) == 1) {
    detection.thresholds <- 10^seq(log10(1e-3), log10(max(data)), length = detection.thresholds)
  }


  if (plot.type == "lineplot") {

    # Calculate the core matrix (prevalence thresholds x abundance thresholds)
    coremat <- core_matrix(x, prevalence.intervals, detection.thresholds)

    res <- core_lineplot(coremat)

  } else if (plot.type == "heatmap") {

    # Here we use taxon x abundance thresholds table indicating prevalences
    data <- otu_table(x)@.Data
    if (taxa_are_rows(x)) {data <- t(data)}
    res <- core_heatmap(data, detection.thresholds = detection.thresholds, palette = palette, min.prevalence = min.prevalence)
    
  }

  p <- res$plot
  p <- p + ggtitle(title)

  if (plot) {
    print(p)
  }

  ret <- list(plot = res$plot, data = res$data,
  	    param = list(prevalence.intervals = prevalence.intervals,
   	    detection.thresholds = detection.thresholds, min.prevalence = min.prevalence))

  ret

}

