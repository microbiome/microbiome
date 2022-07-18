#' @title Global Ecosystem State Variables 
#' @description Global indicators of the ecoystem state, including richness,
#' evenness, diversity, and other indicators
#' @param x A species abundance vector, or matrix (taxa/features x samples)
#' with the absolute count data (no relative abundances), or
#' \code{\link{phyloseq-class}} object
#' @param index Default is ‘NULL’, meaning that all available indices
#' will be included. For specific options, see details.
#' @param zeroes Include zero counts in the diversity estimation.
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
    index <- gsub("richness_", "", index)    
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


#' @title Diversity Index
#' @description Various community diversity indices.
#' @param index Diversity index. See details for options.
#' @param zeroes Include zero counts in the diversity estimation.
#' @inheritParams alpha
#' @return A vector of diversity indices
#' @export
#' @examples
#' data(dietswap)
#' d <- alpha(dietswap, 'shannon')
#' @details
#' By default, returns all diversity indices.
#' The available diversity indices include the following:
#' \itemize{
#' \item{inverse_simpson }{Inverse  Simpson diversity:
#' $1/lambda$ where $lambda=sum(p^2)$ and $p$ are relative abundances.}
#' \item{gini_simpson }{Gini-Simpson diversity $1 - lambda$.
#' This is also called Gibbs–Martin, or Blau index in sociology,
#' psychology and management studies.}
#' \item{shannon }{Shannon diversity ie entropy}
#' \item{fisher }{Fisher alpha; as implemented in the \pkg{vegan} package}
#' \item{coverage }{Number of species needed to cover 50\% of the ecosystem.
#' For other quantiles, apply the function coverage directly.}
#' }
#'   
#' @references
#'
#' Beisel J-N. et al. A Comparative Analysis of Diversity Index Sensitivity.
#' Internal Rev. Hydrobiol. 88(1):3-15, 2003.
#' URL:
#' \url{https://portais.ufg.br/up/202/o/2003-comparative_evennes_index.pdf}
#' 
#' Bulla L. An  index of diversity and its associated diversity measure.
#' Oikos 70:167--171, 1994
#'
#' Magurran AE, McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment
#' (Oxford Univ Press, Oxford), Vol 12.
#'
#' Smith B and Wilson JB. A Consumer's Guide to Diversity Indices.
#' Oikos 76(1):70-82, 1996.
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso dominance, richness, evenness, rarity, alpha
#' @keywords utilities
diversity <- function(x, index="all", zeroes=TRUE) {

    # Only include accepted indices
    index <- tolower(index)
    accepted <- .accepted_diversities()

    # Return all indices
    if (length(index) == 1 && index == "all") {
        index <- accepted
    }
    
    if (!is.null(index)) {
        index <- intersect(index, accepted)
    }
    
    if (!is.null(index) && length(index) == 0) {
        return(NULL)
    }

    tab <- diversities_help(x, index, zeroes)

    if (is.vector(tab)) {
        tab <- as.matrix(tab, ncol=1)
        colnames(tab) <- index
    }

    as.data.frame(tab)

}



diversities_help <- function(x, index="all", zeroes=TRUE) {

    if (length(index)==1 && index=="all") {
        index <- .accepted_diversities()
    }
    
    if (length(index) > 1) {
        tab <- NULL
        nams <- c()
        for (idx in index) {
            newind <- diversities_help(x, index=idx, zeroes=TRUE)
            tab <- cbind(tab, newind)
                if (length(newind)>0) {
                nams <- c(nams, idx)
            }
        }

        colnames(tab) <- nams
        return(as.data.frame(tab))
    }
    
    # Pick data
    otu <- abundances(x)

    if (index == "inverse_simpson") {
        ev <- apply(otu, 2, function(x) {
            inverse_simpson(x)
        })
    } else if (index == "gini_simpson") {
        ev <- apply(otu, 2, function(x) {
            gini_simpson(x)
        })
    } else if (index == "shannon") {

        ev <- apply(otu, 2, function(x) {
            shannon(x)
        })

    } else if (index == "fisher") {
    
        if (length(setdiff(unique(as.vector(otu)%%1), 0)) == 0) {

        estimate <- colSums(otu > 0) > 1
        f <- rep(NA, ncol(otu))
            f[estimate] <- fisher.alpha(otu[,estimate], MARGIN=2)
        ev <- f
        if (!all(estimate)) {
            warning("Some Fisher diversities could not be 
            estimated due to limited observations")
        }

        } else {
    
            warning("Fisher diversity defined only for integers; 
            the OTU table contains non-integers. Fisher not estimated.")
            ev <- NULL
        
        }
    } else if (index == "coverage") {
        ev <- unname(coverage(otu))
    }

    if (!is.null(ev)) {
        names(ev) <- colnames(otu)
    }
    
    ev
    
}


.accepted_diversities <- function () { 
  c("inverse_simpson", "gini_simpson", "shannon",
                    "fisher", "coverage")
}

# x: Species count vector
inverse_simpson <- function(x) {
    
    # Simpson index
    lambda <- simpson_index(x)
    
    # Inverse Simpson diversity
    (1/lambda)
    
}

# x: Species count vector
gini_simpson <- function(x) {
    
    # Simpson index
    lambda <- simpson_index(x)
    
    # Gini-Simpson diversity
    1 - lambda
    
}

simpson_index <- function(x) {
    
    # Relative abundances
    p <- x/sum(x)
    
    # Simpson index
    lambda <- sum(p^2)
    
    lambda
    
}



# x: Species count vector
shannon <- function(x) {

    # Ignore zeroes
    x <- x[x > 0]
    
    # Species richness (number of species)
    S <- length(x)
    
    # Relative abundances
    p <- x/sum(x)
    
    # Shannon index
    (-sum(p * log(p)))
    
}






