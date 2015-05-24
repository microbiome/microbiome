#' Plot abundances with ggplot2
#'
#' @param pseq \code{\link{phyloseq-class}} object
#'
#' @param x A variable to map to the horizontal axis.
#' @param y The variable to map on the vertical axis
#' @param line The variable to map on lines
#' @param color The variable to map on colors
#' @param log10 show y axis on log scale
#' @param title Optional title for the graphic.
#'
#' @return A \code{\link{ggplot}} plot object
#' 
#' @import ggplot2
#' @export
#' @examples # 
#' @keywords utilities

boxplot_abundance <- function (pseq, x, y, line = NULL, color = NULL, log10 = TRUE, title = NULL) {

  xvar <- yvar <- linevar <- colorvar <- NULL

  # Construct example data (df). Ensure that samples are given in same order
  # in metadata and HITChip data.
  # FIXME: this can be removed when official phyloseq package
  # is fixed so as to retain the factor level ordering
  df <- harmonize_fields(sample_data(pseq))
  df$xvar <- factor(as.character(df[[x]]))
  df$yvar <- as.vector(otu_table(pseq)[y, ])

  # Visualize example data with a boxplot
  theme_set(theme_bw(20))
  p <- ggplot(df, aes(x = xvar, y = yvar))
  p <- p + geom_boxplot(fill = "gray") 

  # Add also subjects as lines and points
  if (!is.null(line)) {
    df$linevar <- factor(df[[line]])
    p <- p + geom_line(aes(group = linevar), data = df)
  }

  if (!is.null(color)) {
    df$colorvar <- factor(df[[color]])
    p <- p + geom_point(aes(color = colorvar), data = df, size = 4)
    # Add legend label
    p <- p + guides(color=guide_legend(title=color))
  }
  
  if (log10) {
    p <- p + scale_y_log10()
  }

  # Add axis tests
  p <- p + xlab(x) + ylab("Signal")
  if (is.null(title)) {
    title <- y
  }
  p <- p + ggtitle(title)

  # Plot the image
  p

}
