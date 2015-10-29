#' @title n.phylotypes.per.oligo
#' @description Check number of matching phylotypes for each probe
#' 
#' @param taxonomy oligo - phylotype matching data.frame
#' @param level phylotype level
#' @return number of matching phylotypes for each probe
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
n.phylotypes.per.oligo <- function (taxonomy, level) {
  sapply(split(as.vector(taxonomy[, level]), as.vector(taxonomy[, "oligoID"])), function(x) length(unique(x)))
}


