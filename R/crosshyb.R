#' Cross-hybridization table between multimodal taxa as percentages of shared 
#' probes. The number indicates how many percent of oligos for the row taxon 
#' are also hybridizing the corresponding column taxon.
#'
#' @param tax.level Taxonomic level to investigate
#' @param chip Chip type (e.g. 'HITChip')
#' @param selected.taxa Restrict cross-hyb analysis to the selected groups.
#' @param phylogeny.info phylogeny.info 
#'
#' @return A list containing cross-hybridization table 
#'
#' @examples ch <- CrosshybTable(tax.level = 'L1')
#' @export
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

CrosshybTable <- function(tax.level = "L1", chip = "HITChip", 
    selected.taxa = NULL, 
    phylogeny.info = NULL) {
    
    # Get hylogeny info
    if (is.null(phylogeny.info)) {
        phylogeny.info <- GetPhylogeny(chip, phylogeny.version = "filtered")
    }
    
    # Pick necessary columns
    phi <- phylogeny.info[, c(tax.level, "oligoID")]
    
    # Include only selected groups (if any)
    if (!is.null(selected.taxa)) {
        phi <- phi[phi[[tax.level]] %in% selected.taxa, ]
    }
    
    # Create taxon-oligo mapping matrix
    tax.oligos <- sapply(split(phi, phi[[tax.level]]), function(x) {
        x$oligoID
    })
    tax2oligo <- matrix(0, nrow = length(unique(phi[[tax.level]])), 
                        ncol = length(unique(phi$oligoID)))
    rownames(tax2oligo) <- unique(phi[[tax.level]])
    colnames(tax2oligo) <- unique(phi$oligoID)
    for (tax in names(tax.oligos)) {
        oligos <- tax.oligos[[tax]]
        tax2oligo[tax, oligos] <- 1
    }
    
    # Confusion matrix: how many overlapping oligos between two taxa
    confusion.matrix <- matrix(NA, nrow = length(tax.oligos), 
                               ncol = length(tax.oligos))
    rownames(confusion.matrix) <- names(tax.oligos)
    colnames(confusion.matrix) <- names(tax.oligos)
    for (tax1 in rownames(tax2oligo)) {
        for (tax2 in rownames(tax2oligo)) {
            to <- tax2oligo[c(tax1, tax2), ]
            
            confusion.matrix[tax1, tax2] <- mean(to[tax2, to[tax1, ] == 1])
            confusion.matrix[tax2, tax1] <- mean(to[tax1, to[tax2, ] == 1])
            
        }
    }
    
    confusion.matrix
    
}



#' Cross-hybridization between multimodal taxa as percentages of shared probes. 
#' The number indicates how many percent of oligos for the row taxon are 
#' also hybridizing the corresponding column taxon.
#'
#' @param tax.level Taxonomic level to investigate
#' @param chip Chip type (e.g. 'HITChip')
#' @param selected.taxa Restrict cross-hyb analysis to the selected groups.
#' @param show.plot Produce the plot
#' @param order.rows Order table rows
#' @param order.cols Order table columns
#' @param keep.empty Keep taxa that do not show any cross-hybridization
#' @param rounding Rounding of the cell contents
#' @param phylogeny.info phylogeny.info 
#' @param self.correlations Show self-correlations (always 100 percent); 
#'                          or remove (indicate as 0 percent; default)
#'
#' @return A list containing cross-hybridization table and plot
#'
#' @examples 
#'   # res <- PlotCrosshyb(tax.level = 'L2', rounding = 1, show.plot = FALSE)
#' 
#' @export
#' @import ggplot2
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PlotCrosshyb <- function(tax.level = "L1", chip = "HITChip", 
    selected.taxa = NULL, 
    show.plot = TRUE, order.rows = TRUE, order.cols = TRUE, 
    keep.empty = FALSE, rounding = 1, 
    phylogeny.info = NULL, self.correlations = FALSE) {
    
    # Get crosshyb matrix
    confusion.matrix <- CrosshybTable(tax.level = tax.level, chip = "HITChip", 
        selected.taxa = NULL, 
        phylogeny.info = NULL)
    
    # Remove self-correlations
    if (!self.correlations) {
        diag(confusion.matrix) <- 0
    }
    
    # Focus on selected taxa
    if (!is.null(selected.taxa)) {
        confusion.matrix <- confusion.matrix[selected.taxa, ]
    }
    
    # Remove the taxa that do not have any crosshyb
    if (!keep.empty) {
        confusion.matrix <- confusion.matrix[rowSums(confusion.matrix) > 0, 
            colSums(confusion.matrix) > 0]
    }
    
    # Avoid warnings
    Taxon1 <- Taxon2 <- crosshyb <- NULL
    
    if (length(confusion.matrix) == 0) {
        message(paste("No cross-hybriziation at", tax.level, "level"))
        return(NULL)
    }
    
    # Organize into data frame
    df <- melt(confusion.matrix)
    names(df) <- c("Taxon1", "Taxon2", "crosshyb")
    
    # Switch to percentages
    df[["crosshyb"]] <- 100 * df[["crosshyb"]]
    
    # Order rows and cols
    if (order.rows || order.cols) {
        
        hc <- hclust(as.dist(1 - cor(confusion.matrix)), method = "ward.D")
        colord <- hc$ord
        
        hc <- hclust(as.dist(1 - cor(t(confusion.matrix))), method = "ward.D")
        roword <- hc$ord
        
        if (order.rows) {
            df[["Taxon1"]] <- factor(df[["Taxon1"]], 
                levels = rownames(confusion.matrix)[roword])
        }
        
        if (order.cols) {
            df[["Taxon2"]] <- factor(df[["Taxon2"]], 
                levels = colnames(confusion.matrix)[colord])
        }
    }
    
    # Visualize
    theme_set(theme_bw(15))
    df$labels <- round(df[["crosshyb"]], rounding)
    
    aes <- NULL
    p <- ggplot(df, aes(Taxon2, Taxon1, group = Taxon1))
    p <- p + ggplot2::geom_tile(aes(fill = crosshyb))
    p <- p + ggplot2::geom_text(aes(fill = crosshyb, label = labels, size = 4))
    p <- p + scale_fill_gradient(low = "white", high = "red")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    p <- p + xlab("") + ylab("")
    p <- p + theme(legend.position = "none")
    
    if (show.plot) {
        print(p)
    }
    
    list(data = df[, 1:3], plot = p)
    
}


