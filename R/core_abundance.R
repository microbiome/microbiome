#' @title Core Abundance
#' @description Calculates the community core abundance index.
#' @inheritParams core
#' @return A vector of core abundance indices
#' @export
#' @examples
#' data(dietswap)
#' d <- core_abundance(dietswap, detection=0.1/100, prevalence=50/100)
#' @details The core abundance index gives the relative proportion of the core
#' species (in [0,1]). The core taxa are defined as those that exceed the
#' given population prevalence threshold at the given detection level.
#' @seealso rarity
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_abundance <- function(x, detection=0.1/100, prevalence=50/100) {
    
    # Pick taxa x samples compositional matrix
    xcomp <- abundances(x, transform="compositional")
    
    # Core members
    cm <- core_members(xcomp, detection, prevalence)
    
    # Pick the core and calculate abundance
    if (ncol(xcomp) > 1) {
        colSums(xcomp[cm, ])
    } else {
        sum(xcomp[cm, ])
    }
}



