#' @title Variation line plot
#' @description Plot variation in taxon abundance for many subjects.
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
#'   p <- plot_variation(pseq, "Dialister", tipping.point = 1)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @details Assuming the sample_data(x) has 'subject' field and
#' some subjects have multiple time points.
plot_variation <- function (x, taxon, tipping.point = NULL, lims = NULL, shift = 1e-3) {

  pos <- abundance <- NULL

  m <- sample_data(x)
  otu <- otu_table(x)@.Data
  if (!taxa_are_rows(x)) {otu <- t(otu)}

  d <- otu[taxon, ]

  if (is.null(tipping.point)) {
    tipping.point <- median(d)
  }

  if (is.null(lims)) {
    lims <- round(10*range(d))/10
    lims[[1]] <- lims[[1]] + shift
  }
  
  # Pick subjects with multiple timepoints
  time.subjects <- names(which(table(m$subject) > 1))
  keep <- which(m$subject %in% time.subjects)
  
  ranges <- t(sapply(split(d[keep], as.character(m$subject[keep])), range))
  colnames(ranges) <- c("min", "max")
  
  df <- as.data.frame(ranges)
  df$mid  <- rowMeans(ranges)
  df <- df[order(df$mid),]
  df$pos  <- 1:nrow(df)
  
  # Switches the state
  df$len <- df$max - df$min # Range length
  df$switch <- abs(df$mid - tipping.point) < df$len/2
  dforig <- data.frame(list(abundance = d, subject = m$subject))
  dforig$pos <- df[as.character(dforig$subject), "pos"]

  p <- ggplot()
  p <- p + geom_linerange(data = df, aes(x = pos, ymin = min, ymax = max, color = switch))
  p <- p + scale_color_manual(values = c("black", "red"))
  p <- p + geom_hline(aes(yintercept = tipping.point), linetype = 2, size = 1)
  p <- p + ylab("Abundance")
  p <- p + xlab("Subjects")
  p <- p + guides(color = FALSE)
  #p <- p + ylim(lims[[1]], lims[[2]])
  # TODO add shift also to the axis tick positions to be very exact
  p <- p + ylim(values = range(d) + shift)  
  p <- p + coord_flip()
  p <- p + geom_point(data = dforig, aes(x = pos, y = abundance))
  # p <- p + theme(title = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  p <- p + scale_y_log10(limits = lims)
  p <- p + ggtitle(taxon)
  p
  
}
