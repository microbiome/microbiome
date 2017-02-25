#' @title Core which
#' @description Auxiliary function 
#' @param data phylotypes vs. samples data matrix
#' @param intTr intTr
#' @param prevalenceTr prevalenceTr
#' @return Number of OTUs.
#' @keywords internal
core.which <- function(data, intTr, prevalenceTr) {
    d.bin <- data >= intTr
    prevalences <- rowSums(d.bin)
    nOTUs <- as.numeric(prevalences >= prevalenceTr)
    nOTUs
}


