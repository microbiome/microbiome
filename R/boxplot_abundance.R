#' @title Abundance Boxplot
#' @description Plot phyloseq abundances.
#' @param d \code{\link{phyloseq-class}} object
#' @param x Metadata variable to map to the horizontal axis.
#' @param y OTU to map on the vertical axis
#' @param line The variable to map on lines
#' @param violin Use violin version of the boxplot
#' @param na.rm Remove NAs
#' @param show.points Include data points in the figure
#' @details The directionality of change in paired boxplot is indicated by
#' the colors of the connecting lines.
#' @return A \code{\link{ggplot}} plot object
#' @export
#' @examples
#' data(peerj32)
#' p <- boxplot_abundance(peerj32$phyloseq, x='time', y='Akkermansia',
#'    line='subject')
#' @keywords utilities
boxplot_abundance <-
    function(d, x, y, line=NULL, violin=FALSE, na.rm=FALSE, show.points=TRUE) {
    
    change <- xvar <- yvar <- linevar <- colorvar <- NULL
    pseq <- d
    
    otu <- abundances(pseq)
    # otu <- abundances(pseq) if (!taxa_are_rows(pseq)) {otu <- t(otu)}
    
    # Construct example data (df). Ensure that samples are given in same order
    # in metadata and HITChip data.  FIXME: this can be removed when official
    # phyloseq package is fixed so as to retain the factor level ordering
    
    df <- meta(pseq)
    
    df$xvar <- df[[x]]
    if (!is.factor(df[[x]])) {
        df$xvar <- factor(as.character(df$xvar))
    }
    
    if (y %in% taxa(pseq)) {
        df$yvar <- as.vector(unlist(otu[y, ]))
    } else {
        df$yvar <- as.vector(unlist(sample_data(pseq)[, y]))
    }
    
    if (na.rm) {
        df <- subset(df, !is.na(xvar))
        df <- subset(df, !is.na(yvar))
    }
    
    if (nrow(df) == 0) {
        warning("No sufficient data for plotting available. 
            Returning an empty plot.")
        return(ggplot())
    }
    
    # Factorize the group variable
    df$xvar <- factor(df$xvar)
    
    # Visualize example data with a boxplot
    p <- ggplot(df, aes(x=xvar, y=yvar))
    
    if (show.points) {
        p <- p + geom_point(size=2,
        position=position_jitter(width=0.3), alpha=0.5)
    }
    
    # Box or Violin plot ?
    if (!violin) {
        p <- p + geom_boxplot(fill=NA)
    } else {
        p <- p + geom_violin(fill=NA)
    }
    
    # Add also subjects as lines and points
    if (!is.null(line)) {
        df$linevar <- factor(df[[line]])
        
        # Calculate change directionality
        df2 <- suppressWarnings(df %>%
                                arrange(linevar, xvar) %>%
                                group_by(linevar) %>% 
                                summarise(change=diff(yvar)))
        
        # Map back to data
        df$change <- df2$change[match(df$linevar, df2$linevar)]
        
        # Only show the sign of change for clarity
        df$change <- sign(df$change)
        p <- p + geom_line(data=df,
        aes(group=linevar, color=change), size=1) +
        scale_colour_gradient2(low="blue", mid="black", high="red", 
            midpoint=0, na.value="grey50", guide="none")
    }
    #if (!is.null(color)) {   
    #    df$colorvar <- factor(df[[color]])
    #    # Add legend
    #    # p <- p + geom_point(data=df, aes(color=colorvar), size=4)
    #    # label p <- p + guides(color=guide_legend(title=color))
    #}
    
    # Add axis tests
    p <- p + xlab(x) + ylab(y)
    
    # Return ggplot object
    p
    
}
