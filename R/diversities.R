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
diversities <- function(x, index="all", zeroes=TRUE) {

    .Deprecated("diversity", "The microbiome::diversities function has been 
        replaced by function microbiome::alpha. Update your code accordingly.")
    alpha(x, index="all", zeroes=TRUE)

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
    accepted <- c("inverse_simpson", "gini_simpson", "shannon",
                    "fisher", "coverage")

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

    if (length(index) > 1) {
        tab <- NULL
        for (idx in index) {
            tab <- cbind(tab, diversities_help(x, index=idx, zeroes=TRUE))
        }
        colnames(tab) <- index
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
    
    names(ev) <- colnames(otu)
    
    ev
    
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


