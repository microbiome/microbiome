is.phyloseq <- function (x) {
    length(x) == 1 && is(x) == "phyloseq"
}
