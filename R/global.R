#' @title Global Ecosystem State Variables 
#' @description Global indicators of the ecoystem state, including richness,
#' evenness, diversity, and other indicators
#' @param x A species abundance vector, or matrix (taxa/features x samples)
#' with the absolute count data (no relative abundances), or
#' \code{\link{phyloseq-class}} object
#' @param index Default is ‘NULL’, meaning that all available global indices
#' will be included. For specific options, see details.
#' @param rarity.prevalence Prevalence threshold for determining rare taxa.
#' @param rarity.detection Detection threshold for determining rare taxa.
#' @inheritParams core
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
#' d <- global(dietswap, index='gini')
#' # d <- global(dietswap, index='all')
#'
#' @export
#' @seealso global, dominance, rarity, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
global <- function(x, index="all", rarity.detection = 0.2/100, rarity.prevalence = 20/100) {

    tab <- NULL
    index <- tolower(index)
    index <- gsub("diversities_shannon", "shannon", index)
    index <- unique(index)

    message("Richness")
    if (any(c("all", "richness") %in% index)) {
    
        a <- richness(x)

        if (!is.null(a)) {

            if (is.vector(a)) {
                a <- as.matrix(a, ncol=1)
            }

            colnames(a) <- paste("richness_", colnames(a), sep="")
        
            if (!is.null(tab)) {
                tab <- cbind(tab, a)
            } else {
                tab <- a
            }
        }
    }

    message("Observed (richness 0)")
    if (any(c("all", "observed") %in% index)) {
    
        a <- richness(x, detection = 0)

        if (!is.null(a)) {

            if (is.vector(a)) {
                a <- as.matrix(a, ncol=1)
            }

            colnames(a) <- "observed"
        
            if (!is.null(tab)) {
                tab <- cbind(tab, a)
            } else {
                tab <- a
            }
        }
    }

    message("Diversity")
    a <- diversities(x, index=gsub("diversities_", "", gsub("diversity_", "", index)))
    if (!is.null(a)) {
        if (is.vector(a)) {
            a <- as.matrix(a, ncol=1)
        }
        colnames(a) <- paste("diversities_", colnames(a), sep="")
    
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
    a <- rarity(x, index=gsub("rarity_", "", index), rarity.detection, rarity.prevalence)
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





