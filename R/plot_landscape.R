#' @title Landscape Plot
#' @description Plot abundance landscape ie. sample density in 2D projection landscape
#' @inheritParams get_ordination
#' @param col Variable name to highlight samples (points) with colors
#' @param ... Further arguments to pass
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @details For consistent results, set random seet (set.seed) before function call
#' @examples \dontrun{
#'   data(dietswap)
#'   p <- plot_landscape(dietswap)
#' }
#' @keywords utilities
plot_landscape <- function (x, method = "NMDS", distance = "bray", col = NULL, ...) {

  quiet(proj <- get_ordination(x, method, distance))

  if (is.null(col)) {
    proj$col <- as.factor("black")
  } else {
    proj$col <- proj[[col]]  
  }

  p <- microbiome::densityplot(proj[, 1:2],
                      col = proj$col, legend = T, ...)
  
  p
  
}

