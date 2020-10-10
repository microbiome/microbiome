#' @title Data Transformations for phyloseq Objects
#' @description Standard transforms for \code{\link{phyloseq-class}}.
#' @param x \code{\link{phyloseq-class}} object
#' @param transform Transformation to apply. The options include:
#' 'compositional' (ie relative abundance), 'Z', 'log10', 'log10p',
#' 'hellinger', 'identity', 'clr', or any method from the
#' vegan::decostand function.
#' @param target Apply the transform for 'sample' or 'OTU'.
#' Does not affect the log transform.
#' @param shift A constant indicating how much to shift the baseline
#' abundance (in transform='shift')
#' @param scale Scaling constant for the abundance values when
#' transform = "scale".
#' @return Transformed \code{\link{phyloseq}} object
#' @details In transformation typ, the 'compositional' abundances are returned
#' as relative abundances in [0, 1] (convert to percentages by multiplying
#' with a factor of 100). The Hellinger transform is square root of the
#' relative abundance but instead given at the scale [0,1]. The log10p
#' transformation refers to log10(1 + x). The log10 transformation is applied
#' as log10(1 + x) if the data contains zeroes. CLR transform applies
#' a pseudocount of min(relative abundance)/2 to exact zero relative
#' abundance entries in OTU table before taking logs.
#' @export
#' @examples
#'
#' data(dietswap)
#' x <- dietswap
#'
#' # No transformation
#' xt <- transform(x, 'identity')
#' 
#' # OTU relative abundances
#' # xt <- transform(x, 'compositional')
#' 
#' # Z-transform for OTUs
#' # xt <- transform(x, 'Z', 'OTU')
#'
#' # Z-transform for samples
#' # xt <- transform(x, 'Z', 'sample')
#'
#' # Log10 transform (log10(1+x) if the data contains zeroes)
#' # xt <- transform(x, 'log10')
#'
#' # Log10p transform (log10(1+x) always)
#' # xt <- transform(x, 'log10p')
#'
#' # CLR transform
#' xt <- microbiome::transform(x, 'clr')
#'
#' # Shift the baseline
#' # xt <- transform(x, 'shift', shift=1)
#'
#' # Scale
#' # xt <- transform(x, 'scale', scale=1)
#'
#' @keywords utilities
transform <- function(x, transform = "identity", target = "OTU",
                    shift = 0, scale = 1) {

    y <- NULL
    xorig <- x

    # If x is not a phyloseq object then assume that it is
    # taxa x samples matrix
    if (length(is(x)) == 1 && is.phyloseq(x)) {
        x <- abundances(x)
    }

    if (transform == "relative.abundance") {
        transform <- "compositional"
    }
    
    if (!all(sample(round(prod(dim(abundances(x)))/10)))%%1 == 0) {
        warning("The OTU abundances are not integers. 
        Check that the OTU input data is given as original counts 
        to avoid transformation errors!")
    }
    
    if (transform == "compositional") {
        
        # Assuming taxa x samples matrix
        if (target == "OTU") {
    
            # Minor constant 1e-32 is compared to zero to avoid zero
            # division.  Essentially zero counts will then remain zero
            # and otherwise this wont have any effect.
        
            xt <- apply(x, 2, function(x) {
                x/max(sum(x), 1e-32)
            })
        } else if (target == "sample") {
            xt <- apply(x, 1, function(x) {
                x/max(sum(x), 1e-32)
            })
        }
        
    } else if (transform == "Z") {
        
        # Z transform for sample or OTU
        xt <- ztransform(x, target)
        
    } else if (transform == "clr") {
        
        if (any(abundances(x) < 0)) {
            stop("Non-negative data matrix required for the 
            clr transformation. Exiting.")
        }
        
        # If the data has zeroes, then shift up with a negligible
        # constant to avoid singularities

        xt <- x
        
        # Then transform to compositional data
        xt <- transform(xt, "compositional")
        colnames(xt) <- colnames(x)

        if (any(xt == 0)) {
            v <- as.vector(xt)
            minval <- min(v[v > 0])/2
            xt <- xt + minval
        }
        

        # Pick samples x taxa abundance matrix
        d <- t(apply(xt, 2, function(x) {
            log(x) - mean(log(x))
        }))

        if (nrow(d) == ncol(xt)) {
            rownames(d) <- colnames(xt)
            colnames(d) <- rownames(xt)
        } else {
            colnames(d) <- colnames(xt)
            rownames(d) <- rownames(xt)
        }

        xt <- t(d)

    } else if (transform == "log10") {
        
        # Log transform:
        if (min(x) == 0) {
            warning("OTU table contains zeroes. Using log10(1 + x) transform.")
            # target does not affect the log transform
            xt <- log10(1 + x)
        } else {
            xt <- log10(x)
        }

    } else if (transform == "log10p") {

        xt <- log10(1 + x)
        
    } else if (transform == "identity") {
        
        # No transformation
        xt <- x
        
    } else if (transform == "shift") {
        
        xt <- x + shift
    
    } else if (transform == "scale") {
        
        xt <- scale * x 
        
    } else {
        
        a <- try(xt <- decostand(x, method=transform, MARGIN=2))
        
        if (length(is(a)) == 1 && is(a) == "try-error") {
            xt <- NULL
            stop(paste("Transformation", transform, "not defined."))
        }
    }
    
    xret <- xt
    
    # If the input was phyloseq, then return phyloseq
    if (length(is(xorig)) && is.phyloseq(xorig)) {
        if (taxa_are_rows(xorig)) {
            otu_table(xorig)@.Data <- xret
        } else {
            otu_table(xorig)@.Data <- t(xret)
        }
        xret <- xorig
    }
    
    xret
    
}



#' @title Z Transformation
#' @description Z transform for matrices
#' @details Performs centering (to zero) and scaling (to unit
#'   variance) across samples for each taxa.
#' @param x a matrix
#' @param which margin
#' @return Z-transformed matrix
#' @examples
#' #data(peerj32)
#' #pseqz <- ztransform(abundances(peerj32$phyloseq))
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
ztransform <- function(x, which) {
    
    # Start with log10 transform of the absolute counts
    x <- transform(x, "log10")

    
    if (which == "OTU") {
        
        # Z transform OTUs
        trans <- t(scale(t(x)))
        nullinds <- which(rowMeans(is.na(trans)) == 1)
        
        if (length(nullinds) > 0 & min(x) == 0) {
    
            # warning('Setting undetected OTUs to zero in ztransform')
            # Some OTUs have minimum signal in all samples and scaling
            # gives NA.  In these cases just give 0 signal for these
            # OTUs in all samples
        
            trans[names(nullinds), ] <- 0
        }

        # Use the same matrix format than in original data (taxa x
        # samples or samples x taxa)
    
        xz <- trans

    } else if (which == "sample") {
        
        # Z transform samples
        xz <- apply(x, 2, function(y) {
            (y - mean(y))/sd(y)
        })
        
    }

    xz
    
}




