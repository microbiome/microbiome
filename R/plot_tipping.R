#' @title Variation Line Plot
#' @description Plot variation in taxon abundance for many subjects.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxon Taxonomic group to visualize.
#' @param tipping.point Optional. Indicate critical point for abundance
#' variations to be highlighted.
#' @param lims Optional. Figure X axis limits.
#' @param shift Small constant to avoid problems with zeroes in log10
#' @param xlim Horizontal axis limits
#' @return \code{\link{ggplot}} object
#' @examples
#' data(atlas1006)
#' pseq <- subset_samples(atlas1006, DNA_extraction_method == 'r')
#' pseq <- transform(pseq, 'compositional')
#' p <- plot_tipping(pseq, 'Dialister', tipping.point=1)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @details Assuming the sample_data(x) has 'subject' field and
#' some subjects have multiple time points.
plot_tipping <- function(x, taxon, tipping.point=NULL,
    lims=NULL, shift=0.001, xlim=NULL) {
    
    pos <- abundance <- NULL
    
    m <- meta(x)
    otu <- abundances(x)
    
    d <- otu[taxon, ]

    if (is.null(tipping.point)) {
        message("No tipping point given, indicating the median by dashed line.")
        tipping.point <- median(d)
    }
    
    # Pick subjects with multiple timepoints
    time.subjects <- names(which(table(m$subject) > 1))
    keep <- which(m$subject %in% time.subjects)

    if (length(keep) == 0) {
        warning("No time series in the data.");
        return(ggplot())
    }

    ranges <- t(vapply(split(d[keep], as.character(m$subject[keep])),
        range, c(1,1)))
    colnames(ranges) <- c("min", "max")
    
    df <- as.data.frame(ranges)
    df$mid <- rowMeans(ranges)
    df <- df[order(df$mid), ]
    df$pos <- seq_len(nrow(df))
    
    # Switches the state
    df$len <- df$max - df$min  # Range length
    
    df$switch <- abs(df$mid - tipping.point) < df$len/2
    dforig <- data.frame(list(abundance=d, subject=m$subject))
    dforig$pos <- df[as.character(dforig$subject), "pos"]
    dforig <- subset(dforig, !is.na(pos))

    p <- ggplot()
    p <- p + geom_linerange(data=df,
        aes(x=pos, ymin=min, ymax=max, color=switch))
    p <- p + scale_color_manual(values=c("black", "red"))
    p <- p + geom_hline(aes(yintercept=tipping.point), linetype=2, size=1)
    p <- p + ylab("Abundance")
    p <- p + xlab("Subjects")
    p <- p + guides(color=FALSE)
    p <- p + coord_flip()
    p <- p + geom_point(data=dforig, aes(x=pos, y=abundance))
    p <- p + theme(axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
    
    # Assuming relative abundances
    breaks <- 10^seq(-3, 2, 1)
    if (!is.null(xlim)) {
        breaks <- breaks[breaks < max(xlim)]
    }
    names(breaks) <- as.character(breaks)
    
    if (is.null(xlim)) {
        lims <- c(min(df$min) - 1e-2 * min(df$min),
                max(df$max) + 1e-2 * max(df$max))
        p <- p + scale_y_log10(breaks=breaks,
            labels=names(breaks), limits=lims)
    } else {
        xlim <- sort(xlim)
    xlim[[1]] <- min(min(df$min) - 1e-2 * min(df$min), xlim[[1]])
    xlim[[2]] <- max(max(df$max) + 1e-2 * max(df$max), xlim[[2]])
        p <- p + scale_y_log10(breaks=breaks,
            labels=names(breaks), limits=xlim)
    }
    p <- p + ggtitle(taxon)
    
    p
    
}
