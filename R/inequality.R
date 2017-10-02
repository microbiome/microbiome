#' @title Gini Index
#' @description Calculate Gini indices for a phyloseq object.
#' @inheritParams core
#' @return A vector of Gini indices
#' @examples
#' data(dietswap)
#' d <- inequality(dietswap)
#' @references
#' Relative Distribution Methods in the Social Sciences. Mark S.
#' Handcock and Martina Morris, Springer-Verlag, Inc., New York,
#' 1999. ISBN 0387987789.
#' @seealso diversity, reldist::gini (inspired by that implementation but
#' independently written here to avoid external depedencies)
#' @details Gini index is a common measure for relative inequality in
#' economical income, but can also be used as a community diversity
#' measure. Gini index is between [0,1], and increasing gini index implies
#' increasing inequality.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
inequality <- function(x) {
    
    otu <- abundances(x)
    
    # Gini index for each sample
    do <- apply(otu, 2, function(x) {
        inequality_help(x)
    })
    names(do) <- sample_names(x)
    
    do
    
}



inequality_help <- function(x, w=rep(1, length(x))) {
    # See also reldist::gini for an independent implementation
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    p <- cumsum(w)
    nu <- cumsum(w * x)
    n <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}




