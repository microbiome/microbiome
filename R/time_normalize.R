#' @title Normalize Phyloseq Metadata Time Field
#' @description Shift the time field in phyloseq sample_data such that the
#' first time point of each subject is always 0.
#' @param x phyloseq object. The sample_data(x) should contain the following
#' fields: subject, time
#' @return Phyloseq object with a normalized time field
#' @export
#' @examples
#' data(peerj32)
#' pseq <- time_normalize(peerj32$phyloseq)
time_normalize <- function(x) {
    
    # Store the original input
    xo <- x
    meta <- meta(x)
    
    # Shift the times such that the first time point is always 0
    for (subj in unique(meta$subject)) {
        inds <- which(meta$subject == subj & !is.na(meta$time))
        if (length(inds) > 0) {
            meta[inds, "time"] <- meta[inds, "time"] - min(meta[inds, "time"])
        }
    }
    
    # Fix this to sample metadata in the original input
    sample_data(xo) <- meta
    
    # Return phyloseq with normalized time field
    xo
    
}
