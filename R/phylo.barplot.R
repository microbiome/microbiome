


#' phylo.barplot
#'
#' Description: Barplot for *ITChip sample (across taxa) with higher-level 
#'         taxonomic groups indicated by colours.
#'
#' Arguments:
#'   @param x Data vector across taxa (each element should be named by taxon)
#'   @param color.level Higher-order phylogenetic level to indicate by colors
#'   @param phylogeny.info oligo-phylotype mappings
#'   @param title title
#'   @param plot draw plot TRUE/FALSE
#'   @param sort sort the effects by magnitude
#'
#' Returns:
#'   @return ggplot2 object
#'
#' @export
#' @examples 
#'     data(peerj32)
#'     phylogeny.info <- GetPhylogeny('HITChip', 'filtered'); 
#'     signal <- unlist(peerj32$microbes[1, 1:10]); 
#'     p <- phylo.barplot(signal, 
#'                     color.level = 'L1', 
#'               phylogeny.info = phylogeny.info)
#'                   
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

phylo.barplot <- function(x, color.level = "L1", phylogeny.info = NULL, 
                          title = NULL, plot = TRUE, sort = TRUE) {
    
    if (is.null(phylogeny.info)) {
        warning("phylogeny.info not specified, assuming HITChip phylogeny")
        phylogeny.info <- GetPhylogeny("HITChip", "filtered")
    }
    
    taxa <- names(x)
    tax.lev <- names(which.max(apply(phylogeny.info, 2, function(x) {
        sum(taxa %in% x)
    })))
    
    for (tax.lev in c("oligoID", "species", "L1", "L2")) {
        if (all(taxa %in% phylogeny.info[[tax.lev]])) {
            x.level <- tax.lev
        }
    }
    
    # Collect all into a data.frame
    df <- data.frame(list(taxa = taxa))
    
    # Assign higher-level taxonomic groups
    df[[color.level]] <- unlist(droplevels(levelmap(taxa, level.from = x.level,
                               level.to = color.level, 
        phylogeny.info = phylogeny.info)))
    df[["color.level"]] <- df[[color.level]]
    
    df[["x"]] <- x
    
    # Define colors for L1/L2 groups
    colors <- rainbow(length(unique(df[[color.level]])))
    names(colors) <- as.character(unique(df[[color.level]]))
    
    # Rearrange data.frame
    m <- melt(df)
    
    # Sort by x (ie. change order of factors for plot)
    if (sort) {
        df <- within(df, taxa <- factor(taxa, levels = taxa[order(abs(x))]))
    }
    
    # Plot the image
    p <- ggplot(aes(x = taxa, y = x, fill = color.level), data = df)
    p <- p + scale_fill_manual(
                values = colors[as.character(levels(df[[color.level]]))])
    p <- p + geom_bar(position = "identity", stat = "identity")
    p <- p + theme_bw() + coord_flip()
    p <- p + ylab("Signal") + xlab("") + ggtitle(title)
    p <- p + theme(legend.position = "right")
    p <- p + theme(panel.border = element_rect())
    
    if (plot) {
        print(p)
    }
    
    p
    
} 
