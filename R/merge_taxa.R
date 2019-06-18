#' @title Merge Taxa
#' @description Merge taxonomic groups into a single group.
#' @details In some cases it is necessary to place certain OTUs or other
#' groups into an "other" category. For instance, unclassified groups. This
#' wrapper makes this easy. This function differs from phyloseq::merge_taxa
#' by the last two arguments. Here, in merge_taxa2 the user can specify the
#' name of the new merged group. And the merging can be done based on common
#' pattern in the name.
#' 
#' @param x \code{\link{phyloseq-class}} object
#' @param taxa A vector of taxa names to merge.
#' @param pattern Taxa that match this pattern will be merged.
#' @param name Name of the merged group.
#' @return Modified phyloseq object
#' @examples
#'     data(dietswap)
#'     s <- merge_taxa(dietswap, c())
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
merge_taxa2 <- function(x, taxa = NULL, pattern = NULL, name = "Merged") {

    if (is.null(taxa) && is.null(pattern)) {
        return(x)
    }

    if (!is.null(pattern)) {

        if (!is.null(taxa)) {
            mytaxa <- taxa
        } else {
            mytaxa <- taxa(x)
        }

        if (length(grep(pattern, mytaxa)) == 0) {
            return(x)
        }

        mytaxa <- mytaxa[grep(pattern, mytaxa)]

    } else if (is.null(taxa)) {
        mytaxa <- taxa(x)
    } else {
        mytaxa <- taxa
    }

    x2 <- phyloseq::merge_taxa(x, mytaxa, 1)

    mytaxa <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mytaxa))
    taxa_names(x2) <- gsub(mytaxa[[1]], name, taxa_names(x2))
    tax_table(x2)[nrow(tax_table(x2)),] <- rep(name, ncol(tax_table(x2)))
    
    x2
    
}

