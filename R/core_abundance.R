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
core_abundance <- function(x, detection = 0.1/100, prevalence = 50/100,
    include.lowest = FALSE) {
    
    # Pick taxa x samples compositional matrix
    xcomp <- abundances(transform(x, "compositional"))
    
    # Core members
    cm <- core_members(xcomp, detection, prevalence, include.lowest)

    if (length(cm) == 0) {
        warning("Core is empty with the given abundance and prevalence 
            thresholds. Returning NA for core_abundance. Try to 
            change detection and prevalence parameters.");
        ret <- NA
    }

    # Pick the core and calculate abundance
    if (is.vector(xcomp)) {
        ret <- xcomp # 1 OTU x N samples        
    } else if (length(cm) == 1) {
        ret <- xcomp[cm,]
    } else if (ncol(xcomp) > 1 && length(cm) > 1) {
        ret <- colSums(xcomp[cm, ])        
    } else if (ncol(xcomp) == 1) {
        ret <- sum(xcomp[cm, ])
    }
    ret
}



