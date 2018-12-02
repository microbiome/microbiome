#' @title Collapse Replicate Samples
#' @description Collapse samples, mostly meant for technical replicates.
#' @param x \code{\link{phyloseq-class}} object
#' @param method Collapsing method. Only random sampling ("sample") implemented.
#' @param replicate_id Replicate identifier. A character vector.
#' @param replicate_fields Metadata fields used to determine replicates. 
#' @return Collapsed phyloseq object.
#' @references
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#' data(atlas1006)
#' pseq <- collapse_replicates(atlas1006,
#'         method = "sample",
#'         replicate_fields = c("subject", "time"))
collapse_replicates <- function (x, method = "sample",
    replicate_id = NULL, replicate_fields = NULL) {

    if (!is.null(replicate_fields)) {
        if (!is.null(replicate_id)) {
            stop("Provide only replicate_id 
                OR replicate_fields argument for clarity.")
        }
        replicate_id <- unname(apply(meta(x)[, replicate_fields],
        1, function (x) {paste(x, collapse = "-")}))    
    }

    # Sample names grouped by replicate
    spl <- split(sample_names(x), replicate_id)

    if (method == "sample") {
        # Pick one of the replicates at random
        s <- unname(vapply(spl, function (x) {sample(x, 1)}, "char"))
        x <- prune_samples(s, x)
    }
    # TODO
    # Add averaging of replicates.
    # This is more challenging as we should decide how to merge metadata
    # To be on the safe side, this could be merged only when identical
    # between replicates. Will require checks

    x

}