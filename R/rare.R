#' @title Select Rare Taxa
#' @description Filter the phyloseq object to include only rare taxa.
#' @inheritParams core_members
#' @return Filtered phyloseq object including only rare taxa
#' @references 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   data(peerj32)
#'   rare(peerj32$phyloseq, 200, 20)
rare <- function (x, detection, prevalence) {

  # TODO: add optional renormalization such that the 
  # abundances would sum up to 1 ?

  # Core taxa
  cm <- core_members(x, detection, prevalence)

  # Rare taxa as complement of core taxa
  rt <- setdiff(taxa(x), cm)

  # Pick the subset
  prune_taxa(rt, x)

}




