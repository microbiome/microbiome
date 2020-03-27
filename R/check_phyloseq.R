check_phyloseq <- function (x, fill_na_taxa = TRUE) {

    # Sanity checks for a phyloseq object. Required with some methods.
    if (!taxa_are_rows(x)) {
        x@otu_table <- otu_table(t(otu_table(x)), taxa_are_rows = TRUE)
    }

    if (fill_na_taxa || is.character(fill_na_taxa)) {
        M <- as.matrix(tax_table(x))
        if (!taxa_are_rows(x)) {
            M <- t(M)
        }

        if (!is.character(fill_na_taxa)) {
            fill_na_taxa <- "Unknown"
        }

        M[is.na(M)] <- fill_na_taxa

        x@tax_table <- tax_table(M)

    }

    x

}
