#' @title Global Ecosystem State Variables 
#' @description Global indicators of the ecoystem state, including richness,
#' evenness, diversity, and other indicators
#' @param x A species abundance vector, or matrix (taxa/features x samples)
#' with the absolute count data (no relative abundances), or
#' \code{\link{phyloseq-class}} object
#' @param index Default is ‘NULL’, meaning that all available indices
#' will be included. For specific options, see details.
#' @param zeroes Include zero counts in the diversity estimation.
#' @inheritParams core
#' @return A data.frame of samples x alpha diversity indicators
#' @details This function returns various indices of the ecosystem state.
#' The function is named alpha (global in some previous versions of this
#' package) as these indices can be viewed as measures of
#' alpha diversity. The function uses default choices for detection,
#' prevalence and other parameters for
#' simplicity and standardization. See the individual functions for more
#' options. All indicators from the richness, diversity, evenness,
#' dominance, and rarity functions are available. Some additional measures,
#' such as Chao1 and ACE are available via \code{\link{estimate_richness}}
#' function in the \pkg{phyloseq} package but not included here.
#' The index names are given the prefix richness_, evenness_, diversity_,
#' dominance_, or rarity_ in the output table to avoid confusion between
#' similarly named but different indices (e.g. Simpson diversity and Simpson
#' dominance). All parameters are set to their default. To experiment with
#' different parameterizations, see the more specific index functions
#' (richness, diversity, evenness, dominance, rarity).
#'
#' @examples
#' data(dietswap)
#' d <- alpha(dietswap, index='shannon')
#' # d <- alpha(dietswap, index='all')
#'
#' @export
#' @seealso dominance, rarity, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
alpha <- function(x, index="all", zeroes=TRUE) {

    tab <- NULL
    index <- tolower(index)
    index <- gsub("diversity_shannon", "shannon", index)
    index <- unique(index)

    message("Observed richness")
    if (any(c("all", "observed") %in% index)) {
        a <- richness(x, detection = 0, index = "observed")
        a <- as.matrix(a, ncol=1)
        colnames(a) <- "observed"
        if (!is.null(tab)) {
            tab <- cbind(tab, a)
        } else {
            tab <- a
        }   
    }

    message("Other forms of richness")
    if (any(c("all", "chao1") %in% index)) {
        a <- richness(x, index = "chao1")
        a <- as.matrix(a, ncol=1)
        colnames(a) <- "chao1"

        if (!is.null(tab)) {
            tab <- cbind(tab, a)
        } else {
            tab <- a
        }
    }

    message("Diversity")
    a <- diversity(x, index=gsub("diversity_", "",
        gsub("diversity_", "", index)), zeroes=zeroes)
    if (!is.null(a)) {
        if (is.vector(a)) {
            a <- as.matrix(a, ncol=1)
        }
        colnames(a) <- paste("diversity_", colnames(a), sep="")
    
        if (!is.null(tab)) {
            tab <- cbind(tab, a)
        } else {
            tab <- a
        }

    }

    message("Evenness")
    a <- evenness(x, index=gsub("evenness_", "", index))
    if (!is.null(a)) {
        if (is.vector(a)) {
            a <- as.matrix(a, ncol=1)
        }
        colnames(a) <- paste("evenness_", colnames(a), sep="")

        if (!is.null(tab)) {
            tab <- cbind(tab, a)
        } else {
            tab <- a
        }    

    }

    message("Dominance")
    a <- dominance(x, index=gsub("dominance_", "", index))
    if (!is.null(a)) {
        if (is.vector(a)) {
            a <- as.matrix(a, ncol=1)
        }
        colnames(a) <- paste("dominance_", colnames(a), sep="")

        if (!is.null(tab)) {
            tab <- cbind(tab, a)
        } else {
            tab <- a
        }

    }

    message("Rarity")
    a <- rarity(x, index=gsub("rarity_", "", index))
    if (!is.null(a)) {
        if (is.vector(a)) {
            a <- as.matrix(a, ncol=1)
        }
        colnames(a) <- paste("rarity_", colnames(a), sep="")

        if (!is.null(tab)) {
            tab <- cbind(tab, a)
        } else {
            tab <- a
        }    

    }

    as.data.frame(tab)
}





