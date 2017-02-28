#' @title Data Transformations for phyloseq Objects
#' @description Standard transforms for \code{\link{phyloseq-class}}.
#' @param x \code{\link{phyloseq-class}} object
#' @param transform Transformation to apply. The options include:
#'   'compositional' (ie relative abundance), 'Z', 'log10', 'hellinger',
#'   'identity', 'clr', 'ilr',
#'    or any method from the vegan::decostand function.
#' @param target Apply the transform for 'sample' or 'OTU'.
#'               Does not affect the log transform.
#' @return Transformed \code{\link{phyloseq}} object
#' @details The relative abunance are returned as percentages in [0,
#'   100]. The Hellinger transform is square root of the relative
#'   abundance but instead given at the scale [0,1].
#' @export
#' @examples
#' \dontrun{
#'
#'   # OTU relative abundances
#'   xt <- transform(x, "relative.abundance", "OTU")
#' 
#'   # Z-transform for OTUs
#'   xt <- transform(x, "Z", "OTU")
#'
#'   # Z-transform for samples
#'   xt <- transform(x, "Z", "sample")
#'
#'   # Log10 transform (log(1+x) if the data contains zeroes)
#'   xt <- transform(x, "log10")
#'
#' }
#' @keywords utilities
transform <- function (x, transform = "identity",
                                   target = "OTU") {

  y <- NULL

  if (transform == "relative.abundance") {
    transform <- "compositional"
  }

  if (!all(sample(round(prod(dim(abundances(x)))/10) ))%%1 == 0) {
    warning("The OTU abundances are not integers. 
             Check that the OTU input data is given as original counts 
	     to avoid transform errors!")
  }

  if (transform == "compositional") {
    if (target == "OTU") {
      xt <- transform_sample_counts(x, function (x) {100 * x/sum(x)})
    } else {
      stop(paste("transform not implemented for transform",
      				     transform, "with target", target))
    }
  } else if (transform == "Z") {

    # Z transform for sample or OTU
    xt <- ztransform(x, target)

  } else if (transform == "clr") {

    xt <- x
    if (taxa_are_rows(xt)) {
      a <- t(abundances(transform(xt, "compositional")))
    } else {
      a <- abundances(transform(xt, "compositional"))
    }

    if (!ncol(a) == nsamples(xt)) {stop("Something wrong with clr transform.")}
    
    d <- apply(compositions::clr(a), 2, identity)

    colnames(d) <- sample_names(xt)
    rownames(d) <- taxa(xt)

    xt@otu_table@.Data <- t(d)

  } else if (transform == "log10") {
  
    # Log transform:
    if (min(abundances(x)) == 0) {
      warning("OTU table contains zeroes. Using log10(1 + x) transform.")
      # target does not affect the log transform 
      xt <- transform_sample_counts(x, function(x) log10(1 + x))      
    } else {
      xt <- transform_sample_counts(x, function(x) log10(x))      
    }

  } else if (transform == "identity") {

    # No transform
    xt <- x
    
  } else {
  
    if (target == "OTU") {
    
      xt <- x
      a <- try(xx <- decostand(abundances(xt),
      	                       method = transform, MARGIN = 2))
            
      if (class(a) == "try-error") {
        xt <- NULL
        stop(paste("Transformation", transform, "not defined."))
      }

      if (!taxa_are_rows(xt)) {xx <- t(xx)}
      otu_table(xt)@.Data <- xx

     } else {
    
      stop(paste("Transformation", transform, "not defined for", target))
      
    }    
  }

  xt

}



#' @title Phyloseq Z Transformation
#' @description Z transform phyloseq objects.
#' @details Performs centering (to zero) and scaling (to unit
#'   variance) across samples for each taxa.
#' @param x \code{\link{phyloseq-class}} object 
#' @param which Specify Z transform for "sample" or "OTU"
#' @return Z-transformed phyloseq object
#' @examples \dontrun{
#'   data(peerj32)
#'   pseqz <- ztransform(peerj32$phyloseq)
#' }
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
ztransform <- function (x, which) {

  taxa_are_rows <- y <- NULL

  if (!all(sample(abundances(x), 100)%%1 == 0)) {
    warning("phyloseq object may already have been log transformed - the 
             abundances are not counts - log10 omitted in Z transform. 
	     Perform manually if needed.")
  } else {
    # Start with log10 transform
    x <- transform(x, "log10")
  }
  
  if (which == "OTU") {

    # taxa x samples
    ddd <- abundances(x)

    # Z transform OTUs
    trans <- as.matrix(scale(t(ddd)))

    nullinds <- which(rowMeans(is.na(trans)) == 1)
    if (length(nullinds) > 0 & min(ddd) == 1) {
      warning("Setting undetected OTUs to zero in ztransform")
      # Some OTUs have minimum signal in all samples and scaling gives NA.
      # In these cases just give 0 signal for these OTUs in all samples
      trans[names(which(rowMeans(is.na(trans)) == 1)),] <- 0
    }

    # Use the same matrix format than in original data
    # (taxa x samples or samples x taca)
    xz <- x
    if (taxa_are_rows(x)) { trans = t(trans) }
    otu_table(xz) <- otu_table(trans, taxa_are_rows = taxa_are_rows(x))  

  } else if (which == "sample") {

    # Z transform samples
    xz <- transform_sample_counts(x, function(x) {(y - mean(y))/sd(y) })

  }
  
  xz

}




