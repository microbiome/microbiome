#' @title Visualize Samples of a Microbiota Atlas 
#' @description Show all samples of a microbiota collection, colored by
#' specific factor levels (x axis) and signal (y axis). 
#' @param pseq phyloseq object
#' @param x Sorting variable for X axis and sample coloring
#' @param y Signal variable for Y axis
#' @param ncol Number of legend columns.
#' @return ggplot object
#' @details Arranges the samples based on the given grouping factor (x), and
#' plots the signal (y) on the Y axis. The samples are randomly ordered
#' within each factor level. The factor levels are ordered by standard
#' deviation of the signal (y axis).
#' @references See citation('microbiome');
#' Visualization inspired by Kilpinen et al. 2008,
#' Genome Biology 9:R139. DOI: 10.1186/gb-2008-9-9-r139
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#' data(atlas1006)
#' p <- plot_atlas(atlas1006, 'DNA_extraction_method', 'diversity')
#' p <- plot_atlas(atlas1006, 'DNA_extraction_method', 'Bifidobacterium')
#' @keywords utilities
plot_atlas <- function(pseq, x, y, ncol=2) {
    
    index <- signal <- xvar <- NULL
    
    df <- data.frame(sample_data(pseq)[, x])
    df$xvar <- df[[x]]
    df$xvar <- as.factor(df$xvar)
    
    if (y %in% names(sample_data(pseq))) {
        df$signal <- sample_data(pseq)[[y]]
    } else if (y %in% rownames(abundances(pseq))) {
        df$signal <- as.vector(abundances(pseq)[y, ])
    }
    
    # Randomize sample order to avoid visualization biases
    df <- df[sample(nrow(df)), ]
    
    # Order the factor levels (x variables) by their standard deviation
    df$xvar <- factor(df$xvar,
        levels=names(sort(vapply(split(df$signal, df$xvar), sd, 1))))
    
    # Order by x variable (randomization affects the ordering within
    # each factor level if this is a factor)
    
    df <- df[order(df$xvar), ]
    
    df$index <- seq_len(nrow(df))
    
    p <- ggplot(df, aes(x=index, y=signal, color=xvar))
    p <- p + geom_point()
    p <- p + guides(color=guide_legend(ncol=ncol, title=x))
    p <- p + scale_y_log10()
    p <- p + ggtitle(y)
    p <- p + xlab("Sample index")
    p <- p + ylab(y)
    
    p
    
}

