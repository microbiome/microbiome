#' @title Univariate bimodality plot
#' @description Coloured bimodality plot.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxon Taxonomic group to visualize.
#' @param tipping.point Optional. Indicate critical point for abundance variations to be highlighted.
#' @param lims Optional. Figure X axis limits.
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
#' @details Assuming the sample_data(x) has 'time' field. Uses data only from
#' the baseline (0) time point.
plot_bimodal <- function (x, taxon, tipping.point = NULL, lims = NULL) {

  Abundance <- ..density.. <- ..x.. <- NULL

  m <- sample_data(x)
  otu <- otu_table(x)@.Data
  if (!taxa_are_rows(x)) {otu <- t(otu)}

  d <- otu[taxon, ]

  if (is.null(tipping.point)) {
    tipping.point <- median(d)
  }

  if (is.null(lims)) {
    lims <- round(10*range(na.omit(d)))/10
    lims[[1]] <- lims[[1]] + 1e-3
  }

  lim <- max(abs(d) - tipping.point)

  df <- data.frame(Abundance = d[which(m$time == 0)])

  p <- ggplot(df, aes(x = Abundance, y = ..density.., fill = ..x..)) 
  p <- p + geom_histogram(col = "black", binwidth = .12)
  p <- p + ylab("Frequency") 
  p <- p + xlab("")

  breaks <- c(seq(floor(min(lims)), ceiling(max(lims)), by = 1))
  names(breaks) <- as.character(breaks)

  #p <- p + scale_fill_gradientn(bquote(paste("Signal (", Log[10], ")", sep = "")),
  p <- p + scale_fill_gradientn("Signal",   
                breaks = breaks, 
                colours = c(rep("darkblue", 3), "blue", "white", "red", rep("darkred", 3)), 
                limits = c(-lim + tipping.point, lim + tipping.point),
                labels = names(breaks)
                )

  p <- p + guides(fill = FALSE)
  p <- p + geom_vline(aes(xintercept = tipping.point), linetype = 2, size = 1)
  p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  # p <- p + xlim(lims[[1]], lims[[2]]) 
  p <- p + scale_x_log10(limits = lims) 
  p <- p + ggtitle(taxon) 
  p 

}