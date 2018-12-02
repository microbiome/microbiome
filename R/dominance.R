#' @title Dominance Index
#' @description Calculates the community dominance index.
#' @param index If the index is given, it will override the other parameters.
#' See the details below for description and references of the standard
#' dominance indices. By default, this function returns the Berger-Parker
#' index, ie relative dominance at rank 1.
#' @param rank Optional. The rank of the dominant taxa to consider.
#' @param relative Use relative abundances (default: TRUE)
#' @param aggregate Aggregate (TRUE; default) the top members or not.
#' If aggregate=TRUE, then the sum of relative abundances is returned.
#' Otherwise the relative abundance is returned for the single taxa with
#' the indicated rank. 
#' @inheritParams alpha
#' @return A vector of dominance indices
#' @export
#' @examples
#' data(dietswap)
#' # vector
#' d <- dominance(abundances(dietswap)[,1], rank=1, relative=TRUE)
#' # matrix
#' # d <- dominance(abundances(dietswap), rank=1, relative=TRUE)
#' # Phyloseq object
#' # d <- dominance(dietswap, rank=1, relative=TRUE)
#'
#' @details The dominance index gives the abundance of the most abundant
#' species. This has been used also in microbiomics context
#' (Locey & Lennon (2016)). The following indices are provided:
#' \itemize{
#' \item{'absolute' }{This is the most simple variant, giving the absolute
#' abundance of the most abundant species (Magurran & McGill 2011).
#' By default, this refers to the single most dominant species (rank=1)
#' but it is possible to calculate the absolute dominance with rank n based
#' on the abundances of top-n species by tuning the rank argument.}
#' \item{'relative' }{Relative abundance of the most abundant species.
#' This is with rank=1 by default but can be calculated for other ranks.}
#' \item{'DBP' }{Berger–Parker index, a special case of relative dominance
#' with rank 1; This also equals the inverse of true diversity of the
#' infinite order.}
#' \item{'DMN' }{McNaughton’s dominance. This is the sum of the relative
#' abundance of the two most abundant taxa, or a special case of relative
#' dominance with rank 2}
#' \item{'simpson' }{Simpson's index ($sum(p^2)$) where p are relative
#' abundances has an interpretation as a dominance measure. Also the
#' version ($sum(q * (q-1)) / S(S-1)$) based on absolute abundances q has
#' been proposed by Simpson (1949) but not included here as it is not
#' within [0,1] range, and it is highly correlated with the simpler
#' Simpson dominance. Finally, it is also possible to calculated
#' dominances up to an arbitrary rank by setting the rank argument}
#' \item{'core_abundance' }{Relative proportion of the core species that
#' exceed detection level 0.2\% in over 50\% of the samples}
#' \item{'gini' }{Gini index is calculated with the function
#' \code{inequality}.}
#' }
#'
#' By setting aggregate=FALSE, the abundance for the single n'th most dominant
#' taxa (n=rank) is returned instead the sum of abundances up to that rank
#' (the default).
#'
#' @references
#'
#'   Kenneth J. Locey and Jay T. Lennon.
#'   Scaling laws predict global microbial diversity.
#'   PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#'   Magurran AE, McGill BJ, eds (2011)
#'   Biological Diversity: Frontiers in Measurement and Assessment
#'   (Oxford Univ Press, Oxford), Vol 12
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso coverage, core_abundance, rarity, alpha
#' @keywords utilities
dominance <- function(x, index="all", rank=1, relative=TRUE, aggregate=TRUE) {

    # Only include accepted indices
    if (!is.null(index)) {index <- tolower(index)}
    
    accepted <- tolower(c("DBP", "DMN", "absolute", "relative",
        "simpson", "core_abundance", "gini"))

    # Return all indices
    if (length(index) == 1 && index == "all") {
        index <- accepted
    }

    if (!is.null(index) && !any(index %in% accepted)) {
        return(NULL)
    }

    if (!is.null(index)) {
        index <- intersect(index, accepted)
    }

    tab <- NULL
    if (!is.null(index)) {
        for (idx in index) {
            tab <- cbind(tab,
            dominance_help(x, index=idx, rank=rank,
                relative=relative, 
                aggregate=aggregate))
    }
    } else {
            tab <- cbind(tab,
            dominance_help(x, index=NULL, rank=rank,
                relative=relative, 
                aggregate=aggregate))
    }
    colnames(tab) <- index
    tab <- as.data.frame(tab)
    
    tab

}



dominance_help <- function(x, index="all", rank=1, relative=TRUE,
    aggregate=TRUE) {

    otu <- abundances(x)

    if (is.null(index)) {
        rank <- rank
    } else if (index == "absolute") {
        relative <- FALSE  # Rank=1 by default but can be tuned
    } else if (index %in% c("relative")) {
        relative <- TRUE  # Rank=1 by default but can be tuned
    } else if (index %in% c("dbp")) {
        # Berger-Parker
        rank <- 1
        relative <- TRUE
    } else if (index %in% c("dmn")) {
        # McNaughton's dominance
        rank <- 2
        relative <- TRUE
        aggregate <- TRUE
    } else if (index %in% c("simpson")) {
        tmp <- apply(otu, 2, function(x) {
            simpson_dominance(x)})
        return(tmp)
    } else if (index %in% c("core_abundance")) {
        return(core_abundance(otu, detection=0.2/100, prevalence=50/100))
    } else if (index == "gini") {
        return(inequality(otu))
    }

    if (rank == 1 && relative) {
        index <- "dbp"
    } else if (rank == 2 && relative && aggregate) {
        index <- "dmn"
    }

    if (relative) {
        otu <- apply(otu, 2, function(x) {
            x/sum(x, na.rm=TRUE)
        })
    }
    
    if (!aggregate) {
        do <- apply(otu, 2, function(x) {
            rev(sort(x))[[rank]]
        })
    } else {
        do <- apply(otu, 2, function(x) {
            sum(rev(sort(x))[seq_len(rank)])
        })
    }
    
    names(do) <- colnames(otu)

    if (is.vector(do)) {
        do <- as.matrix(do, ncol=1)
        colnames(do) <- index        
    }

    do
    
}




# x: Species count vector
simpson_dominance <- function(x, zeroes=TRUE) {
    
    if (!zeroes) {
        x[x > 0]
    }
    
    # Relative abundances
    p <- x/sum(x)
    
    # Simpson index (has interpretation as dominance)
    lambda <- sum(p^2)
    
    # More advanced Simpson dominance (Simpson 1949) However let us not use
    # this as it is not in [0,1] and it is very highly correlated with the
    # simpler variant lambda Species richness (number of species)
    # S <- length(x) sum(p * (p - 1)) / (S * (S - 1))

    lambda
    
}

