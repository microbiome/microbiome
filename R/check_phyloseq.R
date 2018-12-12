check_phyloseq <- function (x) {

    # Sanity checks for a phyloseq object. Required with some methods.
    if (!taxa_are_rows(x)) {

        x@otu_table <- otu_table(t(otu_table(x)), taxa_are_rows = TRUE)

    }

    x

}