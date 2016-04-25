#' @title Longitudinal time series plot
#' @description Plot longitudinal time series.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxon Selected taxonomic group.
#' @param subjectID Selected subject.
#' @param tipping.point Optional. Indicate critical point for abundance variations to be highlighted.
#' @param shift Small constant to avoid problems with zeroes in log10
#' @return \code{\link{ggplot}} object
#' @examples 
#'   data("atlas1006")
#'   pseq <- atlas1006
#'   pseq <- transform_phyloseq(pseq, "relative.abundance")
#'   p <- plot_longitudinal(pseq, "Dialister")
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @details Assuming the sample_data(x) has 'subject' field and
#' some subjects have multiple time points.
plot_longitudinal <- function (x, taxon, subjectID = NULL, tipping.point = NULL, shift = 1e-3) {

  m <- sample_data(x)
  
  if (is.null(subjectID) && "subject" %in% colnames(m) && length(unique(m[, "subject"])) == 1) {
    subjectID = as.character(unique(m[, "subject"]))
  }

  if (!is.null(subjectID) && "subject" %in% colnames(m) && length(unique(m[, "subject"])) == 1) {
    x = prune_samples(rownames(m)[which(m[, "subject"] == subjectID)], x)
  } else {
    stop("Specify unique subject ID.")
  }

  x = prune_taxa(taxon, x)

  df = psmelt(x)

  p <- ggplot(data = df, aes(x = time, y = Abundance))
  p <- p + geom_line() + geom_point()
  p <- p + scale_y_log10()
  p <- p + ylab("Abundance")
  p <- p + xlab("Time")
  p <- p + ggtitle(paste(subjectID, taxon, sep = "/"))  

  if (!is.null(tipping.point)) {
    p <- p + geom_hline(aes(yintercept = tipping.point), linetype = 2, size = 1)
  }
  
  p
  
}
