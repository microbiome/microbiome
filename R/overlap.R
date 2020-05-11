#' @title Overlap Measure
#' @description Quantify microbiota 'overlap' between samples.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold.
#' @return Overlap matrix
#' @references
#' Bashan, A., Gibson, T., Friedman, J. et al.
#' Universality of human microbial dynamics.
#' Nature 534, 259–262 (2016). \url{https://doi.org/10.1038/nature18301}
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#' data(atlas1006)
#' o <- overlap(atlas1006, detection = 0.1/100)
overlap <- function(x, detection = 0) {

    # taxa x samples
    a <- abundances(transform(x, "compositional"))
    # samples x samples
    O <- matrix(0, nrow = ncol(a), ncol = ncol(a))
    for (i in seq_len(ncol(a)-1)) {
        for (j in seq(i+1, ncol(a))) {
            O[i, j] <- O[j, i] <- olap(a[,i], a[, j], detection = detection)
        }
    }

    O
    
}


olap <- function (x, y, detection) {

    inds <- which(x > detection & y > detection)
    x <- x[inds]
    y <- y[inds]
    o <- (sum(x) + sum(y))/2

}



