#' @title Global Ecosystem State Variables 
#' @description Global indicators of the ecoystem state, including richness,
#' evenness, diversity, and other indicators
#' @param x A species abundance vector, or matrix (taxa/features x samples)
#' with the absolute count data (no relative abundances), or
#' \code{\link{phyloseq-class}} object
#' @param index Default is ‘NULL’, meaning that all available global indices
#' will be included. For specific options, see details.
#' @inheritParams core
#' @return A data.frame of samples x global indicators
#' @details This function returns global indices of the ecosystem state using
#' default choices for detection, prevalence and other parameters for
#' simplicity and standardization. See the individual functions for more
#' options. All indicators from the richness, diversities, evenness,
#' dominance, and rarity functions are available. Some additional measures,
#' such as Chao1 and ACE are available via \code{\link{estimate_richness}}
#' function in the \pkg{phyloseq} package but not included here.
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
global <- function(x, index="all") {
    
    tab <- NULL
    index <- tolower(index)

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

    a <- diversities(x, index=gsub("diversity_", "", index))
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
    
    if (ncol(tab) > 1) {
        tab <- as.data.frame(tab)
    } else {
        tab <- tab[, 1]
    }

    as.data.frame(tab)
}





