#' @title Plot Cross-Hyb Table
#' @description Cross-hybridization between multimodal taxa as percentages of shared probes. The number indicates how many percent of oligos for the row taxon are also hybridizing the corresponding column taxon.
#' @param tax.level Taxonomic level to investigate
#' @param selected.taxa Restrict cross-hyb analysis to the selected groups.
#' @param show.plot Produce the plot
#' @param order.rows Order table rows
#' @param order.cols Order table columns
#' @param keep.empty Keep taxa that do not show any cross-hybridization
#' @param rounding Rounding of the cell contents
#' @param tax.table tax.table 
#' @param self.correlations Show self-correlations (always 100 percent); 
#'                          or remove (indicate as 0 percent; default)
#' @return A list containing cross-hybridization table and plot
#' @examples 
#'   \dontrun{p <- plot_crosshyb(tax.level = 'L2', rounding = 1, show.plot = FALSE)}
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_crosshyb <- function(tax.level = "L1",
    selected.taxa = NULL, show.plot = TRUE, order.rows = TRUE, order.cols = TRUE, 
    keep.empty = FALSE, rounding = 1, tax.table = NULL, self.correlations = FALSE) {

    ID <- NULL

    # Get crosshyb matrix
    confusion.matrix <- crosshyb_table(tax.level = tax.level, selected.taxa = selected.taxa, tax.table = tax.table)

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
    df <- as.data.frame(confusion.matrix)
    df$ID <- rownames(confusion.matrix)
    df <- gather(df, ID)
    names(df) <- c("Taxon1", "Taxon2", "crosshyb")
    df$crosshyb <- as.numeric(as.character(df$crosshyb))

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
    p <- p + geom_tile(aes(fill = crosshyb))
    p <- p + geom_text(aes(fill = crosshyb, label = labels, size = 4))
    p <- p + scale_fill_gradient(low = "white", high = "red")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    p <- p + xlab("") + ylab("")
    p <- p + theme(legend.position = "none")
    
    if (show.plot) {
      print(p)
    }
    
    list(data = df[, 1:3], plot = p)
    
}


