# Convert data to JSON format for corr_w_scatter visualization
# of correlation matrix linked to scatterplots
# Based on a similar example at 
# http://www.biostat.wisc.edu/~kbroman/D3/corr_w_scatter

#' Description:

#' Arguments:
#'   @param dat data matrix
#'   @param group sample groups (vector or factor)
#'   @param reorder reorder the data
#'
#' Returns:
#'   @return JSON
#'
#' @export
#'
#' @import rjson df2json MASS
#'
#' @examples print("See https://github.com/microbiome/d3/tree/master/corr_w_scatter for an example of the convert4corrwscatter function")
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

convert4corrwscatter <- function(dat, group, reorder=TRUE)
{
  ind <- rownames(dat)
  variables <- colnames(dat)

  if(nrow(dat) != length(group))
    stop("nrow(dat) != length(group)")
  if(!is.null(names(group)) && !all(names(group) == ind))
    stop("names(group) != rownames(dat)")

  if(reorder) {
    ord <- hclust(dist(t(dat)), method="ward")$order
    variables <- variables[ord]
    dat <- dat[,ord]
  }

  # correlation matrix
  corr <- cor(dat, use="pairwise.complete.obs")

  # get rid of names
  dimnames(corr) <- dimnames(dat) <- NULL
  names(group) <- NULL

  # data structure for JSON
  output <- list("ind" = toJSON(ind),
                 "var" = toJSON(variables),
                 "corr" = matrix2json(corr),
                 "dat" =  matrix2json(t(dat)), # columns as rows
                 "group" = toJSON(group))
  paste0("{", paste0("\"", names(output), "\" :", output, collapse=","), "}")
}





