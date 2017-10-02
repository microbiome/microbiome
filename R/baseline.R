#' @title Pick Baseline Timepoint Samples
#' @description Identify and select the baseline timepoint samples in a
#' \code{\link{phyloseq}} object.
#' @param x phyloseq object. Assuming that the sample_data(x) has the fields
#' 'time', 'sample' and 'subject'
#' @param na.omit Logical. Ignore samples with no time point information.
#' If this is FALSE, the first sample for each subject is selected even
#' when there is no time information.
#' @return Phyloseq object with only baseline time point samples selected.
#' @details Arranges the samples by time and picks the first sample for each
#' subject. Compared to simple subsetting at time point zero, this checks
#' NAs and possibility for multiple samples at the baseline, and guarantees
#' that a single sample per subject is selected.
#' @examples
#' data(peerj32)
#' a <- baseline(peerj32$phyloseq)
#' @export
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
baseline <- function(x, na.omit=TRUE) {
    
    # Arrange samples by time, and pick the first sample for each subject
    m <- meta(x)
    
    if (!all(c("time", "sample", "subject") %in% names(m))) {
        stop("The phyloseq sample_data(x) should contain the following fields: 
        time, sample, subject.")
    }
    
    if (na.omit) {
        m <- subset(m, !is.na(time))
    }
    
    suppressWarnings(m <- m %>% arrange(time))
    ss <- m$subject
    s <- m$sample[match(unique(ss), ss)]
    
    xx <- prune_samples(s, x)
    
    xx
}

