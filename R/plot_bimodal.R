#' @title Univariate Bimodality Plot
#' @description Coloured bimodality plot.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxon Taxonomic group to visualize.
#' @param tipping.point Optional. Indicate critical point for abundance variations to be highlighted.
#' @param lims Optional. Figure X axis limits.
#' @param shift Small constant to avoid problems with zeroes in log10
#' @return \code{\link{ggplot}} object
#' @examples 
#'   data("atlas1006")
#'   pseq <- atlas1006
#'   pseq <- subset_samples(pseq, DNA_extraction_method == "r")
#'   pseq <- transform_phyloseq(pseq, "relative.abundance")
#'   p <- plot_bimodal(pseq, "Dialister", tipping.point = 1)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_bimodal <- function (x, taxon, tipping.point = NULL, lims = NULL, shift = 1e-3) {

  Abundance <- ..density.. <- ..x.. <- NULL

  otu <- abundances(x)

  # Add small shift to avoid problems with 0
  # Take log10 to enable useful visualization
  do <- shift + otu[taxon, ]
  do <- log10(do)
  
  d <- do
  
  if (is.null(tipping.point)) {
    tipping.point <- log10(median(10^d) - shift)
  } else {
    tipping.point <- log10(tipping.point - shift)
  }

  if (is.null(lims)) {
    lims <- range(na.omit(d))
  } else {
    lims <- lims
  }
  lims[[1]] <- floor(100 * lims[[1]])/100
  lims[[2]] <- ceiling(100 * lims[[2]])/100
  lims <- lims + shift
  lim <- max(abs(lims)) 

  breaks <- c(seq(floor(min(lims)), ceiling(max(lims)), by = 1))
  names(breaks) <- as.character(10^breaks)

  lims2 <- c(log10(tipping.point) - lim - .2,
             log10(tipping.point) + lim + .2)

  # Data
  # bquote(paste("Signal (", Log[10], ")", sep = ""))  
  df <- data.frame(Abundance = d)
  p <- ggplot(df, aes(x = Abundance, y = ..density.., fill = ..x..)) +
         geom_histogram(col = "black", binwidth = .12) +
	 ylab("Frequency") +
	 xlab("") +
	 scale_fill_gradientn("Signal",   
                breaks = breaks - 10^tipping.point, 
                colours = c(rep("darkblue", 3), "blue", "white", "red", rep("darkred", 3)),
                labels = names(breaks),
		limits = lims2
                ) +
	guides(fill = FALSE) +
	geom_vline(aes(xintercept = tipping.point), linetype = 2, size = 1) +
	theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
	scale_x_continuous(breaks = breaks, labels = names(breaks), limits = lims) +
	ggtitle(taxon)

  p 

}