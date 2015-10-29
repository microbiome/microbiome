#' @title rda_physeq
#' @description RDA for phyloseq objects.
#'           Based on the \code{\link{rda}} function 
#' 	     from the \pkg{vegan} package.
#'
#' @param x \code{\link{phyloseq-class}} object
#' @param varname Variable to apply in RDA visualization.
#' @param scale See help(rda)
#' @param na.action See help(rda)
#'
#'   @return rda result. See help(vegan::rda)
#'
#' @export
#' @importFrom vegan rda
#'
#' @examples #
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rda_physeq <- function (x, varname, scale = FALSE, na.action = na.fail) {

  # Microbiota profiling data (44 samples x 130 bacteria)
  otu <- t(otu_table(x)@.Data)

  # Sample annotations
  annot <- sample_data(x)[[varname]]

  # Run RDA
  rdatest <- rda(otu ~ annot, scale = scale, na.action = na.action) 
  
  rdatest

}


