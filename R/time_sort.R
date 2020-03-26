#' @title Temporal Sorting Within Subjects
#' @description Within each subject, sort samples by time and calculate
#' distance from the baseline point (minimum time).
#' @param x A metadata data.frame including the following columns:
#' time, subject, sample, signal. Or a phyloseq object.
#' @return A list with sorted metadata (data.frame) for each subject.
#' @references See citation('microbiome')
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#' data(atlas1006)
#' pseq <- subset_samples(atlas1006, DNA_extraction_method == "r")
#' ts <- time_sort(meta(pseq))
#' @keywords utilities
time_sort <- function(x) {
    
    if (is.phyloseq(x)) {
        x <- meta(x)
    }
    
    # Keep only samples with time point info
    x <- subset(x, !is.na(time))

    if (nrow(x) == 0) {
        return(NULL)
    }

    if (is.null(x$signal)) {
        x$signal <- rep(NA, nrow(x))
    }
    
    # Pick data for each subject separately
    spl <- split(x, as.character(x$subject))

    # Keep only subjects with multiple time points
    spl <- spl[names(which(vapply(spl, function(s) {
        length(unique(na.omit(s$time)))
    }, 1) > 1))]
    
    # Ignore NA times
    spl <- lapply(spl, function(s) {
        s[!is.na(s$time), ]
    })
    
    tabs <- list()
    cnt <- 0
    for (subj in names(spl)) {
        
        times <- as.numeric(spl[[subj]]$time)
        signal <- as.numeric(spl[[subj]]$signal)
        mintime <- which.min(times)
        
        # Shift in time from first time point
        spl[[subj]]$time <- (times - times[[mintime]])
        
        # Shift in signal from first time point
        spl[[subj]]$shift <- (signal - signal[[mintime]])
        
        # Store
        tabs[[subj]] <- spl[[subj]]
    }
    
    tabs
    
}




