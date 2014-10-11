#' Description: Calculate distance matrix between the _columns_ of the 
#' input matrix. Can prduce correlation-based distance matrices, otherwise
#' uses the standard 'dist' function.
#'
#' Arguments:
#'   @param x data matrix
#'   @param method distance method
#'   @param ... other arguments to be passed
#'
#' Returns:
#'   @return distance object
#'
#' @export
#'
#' @examples 
#'   data(peerj32)
#'   d <- distance.matrix(peerj32$microbes[1:10, 1:3])
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

distance.matrix <- function(x, method = "pearson", ...) {
    if (method %in% c("pearson", "spearman")) {
        cmat <- as.dist((1 - cor(x, use = "pairwise.complete", 
               method = method)))
    } else {
        cmat <- dist(x, method = method, ...)
    }
    cmat
}


