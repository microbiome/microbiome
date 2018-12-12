check_phyloseq <- function (x, fill_na_taxa = TRUE) {

    # Sanity checks for a phyloseq object. Required with some methods.
    if (!taxa_are_rows(x)) {
        x@otu_table <- otu_table(t(otu_table(x)), taxa_are_rows = TRUE)
    }

    if (fill_na_taxa) {
      M <- as.matrix(tax_table(x))
      M[is.na(M)] <- "Unknown"
      x@tax_table <- tax_table(M)
    }
    
    x

}




#aggregate_na_level <- function (x) {
#  tx <- rownames(tax_table(x))[which(is.na(tax_table(x)[, level]))]
#  x <- remove_taxa(taxa, x)
#}