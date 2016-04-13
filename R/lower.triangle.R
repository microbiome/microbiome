#' @title Lower triangle of a matrix
#' @description Get lower triangle of a square matrix.
#' as a numeric vector such that
#' row-by-row, picking elements in the order
#' 2,1;3,1;3,2;4,1,...
#'        
#' @param mat data matrix
#'
#' @return lower triangle as vector 
#'
#' @export
#' @examples 
#'   mat <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
#'   vec <- lower.triangle(mat)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
lower.triangle <- function(mat) {
    
  # TODO is this easy replace with standard R functions ?

    elements <- c()
    nr <- dim(mat)[[1]]
    nc <- dim(mat)[[2]]
    
    for (i in 2:nr) {
        for (j in 1:(i - 1)) {
            elements <- c(elements, mat[i, j])
        }
    }
    elements
} 
