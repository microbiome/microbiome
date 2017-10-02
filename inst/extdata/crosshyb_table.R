#' @title Cross Hybridization Table
#' @description Cross-hybridization table between multimodal taxa as
#'              percentages of shared probes.
#' The number indicates how many percent of oligos for the row taxon are
#' also hybridizing the corresponding column taxon.
#' @param tax.level Taxonomic level to investigate
#' @param selected.taxa Restrict cross-hyb analysis to the selected groups.
#' @param tax.table tax.table 
#' @return A list containing cross-hybridization table 
#' @examples \dontrun{ch <- crosshyb_table(tax.level = 'L1')}
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
crosshyb_table <- function(tax.level = "L1", selected.taxa = NULL, tax.table) {

    # Pick necessary columns
    phi <- tax.table[, c(tax.level, "oligoID")]

    # Include only selected groups (if any)
    if (!is.null(selected.taxa)) {
        phi <- phi[phi[, tax.level] %in% selected.taxa, ]
    }

    # Create taxon-oligo mapping matrix
    spl <- split(as.character(phi[, "oligoID"]), as.character(phi[, tax.level]))

    tax.oligos <- sapply(spl, function(x) {x})

    tax2oligo <- matrix(0, nrow = length(unique(phi[, tax.level])), 
                        ncol = length(unique(phi[, "oligoID"])))
    rownames(tax2oligo) <- unique(phi[, tax.level])
    colnames(tax2oligo) <- unique(phi[, "oligoID"])
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


