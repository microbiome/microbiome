#' @title Visualize OTU Core
#' @description Core visualization (2D).
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range,
#'          or a scalar indicating the number of intervals in the data range.
#' @param plot.type Plot type ('lineplot' or 'heatmap')
#' @param colours colours for the heatmap
#' @param min.prevalence If minimum prevalence is set, then filter out those
#'    rows (taxa) and columns (detection thresholds) that never exceed this
#'    prevalence threshold. This helps to zoom in on the actual core region
#'    of the heatmap. Only affects the plot.type = 'heatmap'.
#' @param taxa.order Ordering of the taxa.
#' @param horizontal Logical. Horizontal figure.
#' @return A list with three elements: the ggplot object and the data.
#'         The data has a different form for the lineplot and heatmap.
#'         Finally, the applied parameters are returned.
#' @examples 
#'   data(atlas1006)
#'   pseq <- atlas1006
#'   p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10),
#'                        detection.thresholds = c(0, 10^(0:4)))
#' @export 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_core <- function(x, 
		   prevalence.intervals = seq(5, 100, 5), 		   
		   detection.thresholds = 20,
		   plot.type = "lineplot",
		   colours = gray(seq(0,1,length=5)), 		   
		   min.prevalence = NULL,
		   taxa.order = NULL,
		   horizontal = FALSE) {

  if (length(detection.thresholds) == 1) {
    detection.thresholds <- 10^seq(log10(1e-3),
      log10(max(abundances(x), na.rm = T)),
      length = detection.thresholds)
  }

  if (plot.type == "lineplot") {

    # Calculate the core matrix (prevalence thresholds x abundance thresholds)
    coremat <- core_matrix(x, prevalence.intervals, detection.thresholds)

    res <- core_lineplot(coremat)

  } else if (plot.type == "heatmap") {

    # Here we use taxon x abundance thresholds table indicating prevalences
    res <- core_heatmap(abundances(x),
    	                detection.thresholds = detection.thresholds,
			colours = colours, min.prevalence = min.prevalence,
			taxa.order = taxa.order)
    
  }

  p <- res$plot + ggtitle("Core")

  if (horizontal) {
    p <- p + coord_flip() + theme(axis.text.x = element_text(angle = 90))
  }

  #ret <- list(plot = p, data = res$data,
  # 	    param = list(prevalence.intervals = prevalence.intervals,
  # 	    detection.thresholds = detection.thresholds,
  #	    min.prevalence = min.prevalence))

  p

}

