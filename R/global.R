#' @title Global Ecosystem State Variables 
#' @description Global indicators of the ecoystem state, including richness,
#' evenness, diversity, and other indicators
#' @inheritParams alpha
#' @return A data.frame of samples x global indicators
#' @details This function returns global indices of the ecosystem state using
#' default choices for detection, prevalence and other parameters for
#' simplicity and standardization. See the individual functions for more
#' options. All indicators from the richness, diversities, evenness,
#' dominance, and rarity functions are available. Some additional measures,
#' such as Chao1 and ACE are available via \code{\link{estimate_richness}}
#' function in the \pkg{phyloseq} package but not included here.
#' The index names are given the prefix richness_, evenness_, diversities_,
#' dominance_, or rarity_ in the output table to avoid confusion between
#' similarly named but different indices
#' (e.g. Simpson diversity and Simpson dominance).
#'
#' @examples
#' data(dietswap)
#' d <- alpha(dietswap, index='observed')
#'
#' @export
#' @seealso global, dominance, rarity, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
global <- function(x, index="all") {

    .Deprecated("alpha", msg = "The microbiome::global function is deprecated. 
        Use the function microbiome::alpha instead.")
    alpha(x, index) 

}





