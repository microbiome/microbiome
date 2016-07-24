#' @title Abundance Boxplot
#' @description Plot phyloseq abundances.
#' @param pseq \code{\link{phyloseq-class}} object
#' @param x Metadata variable to map to the horizontal axis.
#' @param y OTU to map on the vertical axis
#' @param line The variable to map on lines
#' @param color The variable to map on colors
#' @param log10 show y axis on log scale
#' @details The directionality of change in paired boxplot is indicated by the colors of the connecting lines.
#' @return A \code{\link{ggplot}} plot object
#' @export
#' @examples
#'   data(peerj32)
#'   p <- boxplot_abundance(peerj32$phyloseq, x = "time", y = "Akkermansia",
#'       	  		    line = "subject", color = "gender")
#' @keywords utilities
boxplot_abundance <- function (pseq, x, y, line = NULL, color = NULL, log10 = TRUE) {

  change <- xvar <- yvar <- linevar <- colorvar <- NULL

  # Construct example data (df). Ensure that samples are given in same order
  # in metadata and HITChip data.
  # FIXME: this can be removed when official phyloseq package
  # is fixed so as to retain the factor level ordering
  #df <- harmonize_fields(sample_data(pseq))
  df <- sample_data(pseq)
  df$xvar <- factor(as.character(df[[x]]))
  df$yvar <- as.vector(otu_table(pseq)[y, ])

  # Visualize example data with a boxplot
  theme_set(theme_bw(20))
  p <- ggplot(df, aes(x = xvar, y = yvar))
  p <- p + geom_boxplot(fill = "gray") 

  # Add also subjects as lines and points
  if (!is.null(line)) {
    df$linevar <- factor(df[[line]])

    # Calculate change directionality
    df2 <- suppressWarnings(df %>% arrange(linevar, xvar) %>%
    	   		           group_by(linevar) %>%
				   summarise(change = diff(yvar)))
    
    # Map back to data
    df$change <- df2$change[match(df$linevar, df2$linevar)]
    # Log10 for line colors
    #df$change <- sign(df$change) * log10(1 + abs(df$change))
    # Only show the sign of change for clarity
    df$change <- sign(df$change) 
    p <- p + geom_line(data = df, aes(group = linevar, color = change), size = 1)

    p <- p + scale_colour_gradient2(low = "blue",
      	     			    mid = "black",
				    high = "red",
      	                            midpoint = 0,
				    na.value = "grey50",
				    guide = "none")
   
  }


  if (!is.null(color)) {
  
    df$colorvar <- factor(df[[color]])
    # p <- p + geom_point(data = df, aes(color = colorvar), size = 4)

    # Add legend label
    # p <- p + guides(color = guide_legend(title = color))
    
  }


  if (log10) {
    p <- p + scale_y_log10()
  }

  # Add axis tests
  p <- p + xlab(x) + ylab("Signal")
  if (is.null(title)) {
    title <- y
  }

  # Return ggplot object
  p

}
