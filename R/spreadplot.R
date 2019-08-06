#' @title Abundance Spread Plot
#' @description Visualize abundance spread for OTUs
#' @param x \code{\link{phyloseq-class}} object; or a data.frame with fields
#'        "otu" (otu name); "sample" (sample name); and "abundance"
#'        (otu abundance in the given sample)
#' @param trunc Truncate abundances lower than this to zero
#' @param alpha Alpha level for point transparency
#' @param width Width for point spread
#' @return ggplot2 object
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @examples
#' data(dietswap)
#' p <- spreadplot(transform(dietswap, "compositional"))
#' @keywords utilities
spreadplot <- function (x, trunc = 0.001/100, alpha = 0.15, width = 0.35) {

    otu <- sample <- abundance <- NULL

    # df: data.frame with fields "otu" (otu name);
    # "sample" (sample name); and "abundance"
    # (otu abundance in the given sample)

    df <- melt(abundances(x)) 
    names(df) <- gsub("Var1", "otu", names(df))
    names(df) <- gsub("Var2", "sample", names(df))
    names(df) <- gsub("value", "abundance", names(df))

    o <- df %>% group_by(otu) %>%
            summarise(median = median(abundance),
                mean = mean(abundance)) %>%
        arrange(median) %>%
        mutate(otu = factor(otu, unique(otu)))

    top <- rev(as.character(levels(o$otu))) #[1:50]
    df <- subset(df, otu %in% top)
    df$otu <- factor(df$otu, levels=rev(top))
    brs <- c(10^(-rev(seq(0:3))), 0.5)
    # Truncate abundances below the threshold
    df$abundance[df$abundance < trunc] <- trunc

    p <- ggplot(df, aes(x=otu, y=abundance)) +
        # geom_boxplot(fill = "gray") +
        geom_jitter(alpha = alpha, width = width) +       
        scale_y_continuous(# labels=scales::percent,
                        labels = paste0(100 * brs, "%"),        
                        trans  = "log10",
                        breaks = brs
            #limits = c(1e-5, 1)
            ) +
        coord_flip()

    p

}