#' @title Filter Prevalent Taxa
#' @description Filter the phyloseq object to include only prevalent taxa.
#' @inheritParams core_members
#' @return Filtered phyloseq object including only prevalent taxa
#' @references 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   data(peerj32)
#'   core(peerj32$phyloseq, 200, 20)
core <- function (x, detection, prevalence) {
  # Was: filter_prevalent

  # TODO: add optional renormalization such that the core member
  # abundances would sum up to 100 ?

  taxa <- core_members(x, detection, prevalence)
  prune_taxa(taxa, x)
}




