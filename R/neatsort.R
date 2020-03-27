#' @title Neatmap Sorting
#' @description Sort samples or features based on the neatmap approach. 
#' @param x \code{\link{phyloseq-class}} object or a matrix
#' @param target For \code{\link{phyloseq-class}} input, the target is either
#' 'sites' (samples) or 'species' (features) (taxa/OTUs); for matrices,
#' the target is 'rows' or 'cols'.
#' @param method Ordination method. See \code{\link{ordinate}}
#' from \pkg{phyloseq} package.
#' For matrices, only the NMDS method is available.
#' @param distance Distance method. See \code{\link{ordinate}}
#' from \pkg{phyloseq} package.
#' @param first Optionally provide the name of the first sample/taxon to
#' start the ordering (the ordering is cyclic so we can start at any
#' point). The choice of the first sample may somewhat affect the
#' overall ordering.
#' @param ... Arguments to be passed.
#' @return Vector of ordered elements
#' @export
#' @examples
#' data(peerj32)
#' pseq <- peerj32$phyloseq
#' # For Phyloseq
#' sort.otu <- neatsort(pseq, target='species')
#' # For matrix
#' # sort.rows <- neatsort(abundances(pseq), target='rows')
#'
#' @references This function is partially based on code derived from the
#' \pkg{phyloseq} package. For the original
#' neatmap approach for heatmap sorting, see (and cite):
#' Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map
#' alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' @details This function borrows elements from the heatmap implementation in
#' the \pkg{phyloseq} package. The row/column sorting is there not available
#' as a separate function at present, however, hindering reuse in other tools.
#' Implemented in the microbiome package to provide an independent method for
#' easy sample/taxon reordering for phyloseq objects.
#' @keywords utilities
neatsort <- function(x, target, method="NMDS", distance="bray",
    first=NULL,  ...) {
    
    # TODO harmonize completely for matrices vs phyloseqs
    xo <- x
    
    # Remame matrix targets as in phyloseq
    target <- gsub("columns", "sites", target)
    target <- gsub("cols", "sites", target)
    target <- gsub("rows", "species", target)
    target <- gsub("samples", "sites", target)
    target <- gsub("features", "species", target)

    if (length(is(x)) == 1 && is(x) %in% c("phyloseq", "otu_table")) {

        samples <- sample_names(x)
        features <- taxa(x)
        
        # Capture the output to keep the screen clean
        junk <- capture.output(
            ord <- ordinate(x, method = method, distance = distance, ...),
                file=NULL)
        
    } else {

        samples <- colnames(x)
        features <- rownames(x)
        
        if (!method == "NMDS") {
            warning("Only NMDS dissimilarity implemented for matrices. 
                    Using NMDS for sorting.")
            method <- "NMDS"
        }
        
        if (target == "sites") {
            x <- t(x)
        }
        
        # Neatmap sorting for matrices with NMDS Order Capture the output
        # to keep the screen clean
        if (distance %in% c("euclidean")) {
            d <- dist(x, distance)
        } else {
            d <- vegdist(x, distance)
        }
        
        junk <- capture.output(ord <- metaMDS(d, wascores=FALSE,
        autotransform=FALSE, noshare=FALSE), file=NULL)
    }
    
    # ------------------------------------------------------------------

    # Order items with the NeatMap procedure Reorder by the angle in radial
    # coordinates on the 2-axis plane.
    DF <- NULL
    
    # Define new sample ordering based on the ordination Quick fix: the scores
    # function fails otherwise.
    disp.target <- target
    if (target == "species" && !is.phyloseq(xo)) {
        disp.target <- "sites"
    }
    
    tmp <- try({
        DF <- scores(ord, choices=c(1, 2), display=disp.target)
    }, silent=TRUE)
    
    if (inherits(tmp, "try-error")) {
        warning(paste("Order failed with ", target, ". 
            Using default ordering.", 
            sep=""))
    }
    
    if (!is.null(DF)) {
        # If the score accession worked, replace order
        if (target == "sites") {
            ordering <- samples[order(radial_theta(DF))]
        } else if (target == "species") {
            ordering <- features[order(radial_theta(DF))]
        } else {
            stop("Target should be either sites, cols, species or rows")
        }
    } else if (length(DF) > 1) {
        
        if (target == "sites") {
            ordering <- samples[order(DF)]  # 1:nsamples(x)
        } else if (target == "species") {
            ordering <- features[order(DF)]  # 1:ntaxa(x)
        } else {
            stop("Target should be either sites or species")
        }
    } else {
        if (target == "sites") {
            ordering <- samples
        } else if (target == "species") {
            ordering <- features
        } else {
            stop("Target should be either sites or species")
        }
    }
    
    # Determine the starting item (OTU or sample)
    if (!is.null(first)) {
        ordering <- chunk_reorder(ordering, first)
    }
    
    ordering
    
}



#' @title Radial Theta Function
#' @description Adapted from \pkg{NeatMap} and \pkg{phyloseq} packages
#'   but not exported and hence not available via phyloseq. Completely
#'   rewritten to avoid license conflicts. Vectorized to gain efficiency;
#'   only calculates theta and omits r.
#' @param x position parameter
#' @return theta
#' @keywords internal
radial_theta <- function(x) {
    
    x <- as(x, "matrix")
    theta <- atan2((x[, 2] - mean(x[, 2])), (x[, 1] - mean(x[, 1])))
    names(theta) <- rownames(x)
    
    theta
    
}

#' @title Chunk Reorder
#' @description Chunk re-order a vector so that specified newstart is first.
#' Different than relevel.
#' @keywords internal
#' @details Borrowed from \pkg{phyloseq} package as needed here and not
#' exported there. Rewritten.
#' @return Reordered x
#' @examples 
#' # Typical use-case
#' # chunk_reorder(1:10, 5)
#' # # Default is to not modify the vector
#' # chunk_reorder(1:10)
#' # # Another example not starting at 1
#' # chunk_reorder(10:25, 22)
#' # # Should silently ignore the second element of `newstart`
#' # chunk_reorder(10:25, c(22, 11))
#' # # Should be able to handle `newstart` being the first argument already
#' # # without duplicating the first element at the end of `x`
#' # chunk_reorder(10:25, 10)
#' # all(chunk_reorder(10:25, 10) == 10:25)
#' # # This is also the default
#' # all(chunk_reorder(10:25) == 10:25)
#' # # An example with characters
#' # chunk_reorder(LETTERS, 'G') 
#' # chunk_reorder(LETTERS, 'B') 
#' # chunk_reorder(LETTERS, 'Z') 
#' # # What about when `newstart` is not in `x`? Return x as-is, throw warning.
#' # chunk_reorder(LETTERS, 'g') 
chunk_reorder <- function(x, newstart=x[[1]]) {
    
    pivot <- match(newstart[1], x, nomatch=NA)
    
    # If pivot is NA, then warn and return x as-is
    if (is.na(pivot)) {
        warning("`newstart` argument not found from `x`. 
        - returning `x` with no reordering.")
        newx <- x
    } else {
        newx <- c(tail(x, {
            length(x) - pivot + 1
        }), head(x, pivot - 1L))
    }
    
    newx
    
}
