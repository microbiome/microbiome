#' @title Univariate bimodality plot
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

  otu <- taxa_abundances(x)

  # Add small shift to avoid problems with 0
  # Take log10 to enable useful visualization
  do <- log10(shift + otu[taxon, ])
  d <- do
  
  if (is.null(tipping.point)) {
    tipping.point <- median(10^d) - shift
  } else {
    tipping.point <- tipping.point - shift
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
  df <- data.frame(Abundance = d)
  p <- ggplot(df, aes(x = Abundance, y = ..density.., fill = ..x..)) 
  p <- p + geom_histogram(col = "black", binwidth = .12)
  p <- p + ylab("Frequency") 
  p <- p + xlab("")
  
  # bquote(paste("Signal (", Log[10], ")", sep = ""))
  p <- p + scale_fill_gradientn("Signal",   
                breaks = breaks - tipping.point, 
                colours = c(rep("darkblue", 3), "blue", "white", "red", rep("darkred", 3)),
                labels = names(breaks),
		limits = lims2
                )

  p <- p + guides(fill = FALSE)
  p <- p + geom_vline(aes(xintercept = log10(tipping.point)), linetype = 2, size = 1)
  p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p <- p + scale_x_continuous(breaks = breaks, labels = names(breaks), limits = lims)   
  p <- p + ggtitle(taxon)

  p 

}