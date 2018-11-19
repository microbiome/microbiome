#' @title Time Split
#' @description For each subject, return temporally consecutive sample pairs
#'   together with the corresponding time difference.
#' @param x \pkg{phyloseq} object.
#' @return data.frame
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples
#' data(atlas1006)
#' x <- timesplit(subset_samples(atlas1006,
#'   DNA_extraction_method == 'r' & sex == "male"))
#' @keywords utilities
timesplit <- function(x) {
    
    x0 <- x
    df <- sample_data(x0)

    # Ensure time is numeric
    df$time <- as.numeric(as.character(df$time))
    
    # Remove subjects with only one measurement
    df <- df[df$subject %in% names(which(table(df$subject) > 1)), ]
    
    # Split data by subject
    spl <- split(df, as.character(df$subject))
    spl.list <- list()
    dfis <- NULL    
    for (i in seq_len(length(spl))) {
        # Ensure the measurements are ordered in time
        spl.list[[i]] <- list(spl=spl[[i]][order(spl[[i]]$time), ],
        time.difs=diff(spl[[i]]$time))
        
        spli <- spl.list[[i]]$spl
        time.difs <- spl.list[[i]]$time.difs
    n <- nrow(spl.list[[i]]$spl)
    subject <- unique(spli$subject)
    sample1 <- spl.list[[i]]$spl$sample[seq_len(n-1)]
    sample2 <- spl.list[[i]]$spl$sample[seq(2, n)]    

        # Organize into data frame
        dfi <- data.frame(
        subject = rep(subject, n-1),
        sample1 = sample1,
        sample2 = sample2,
        time=time.difs)
        
        # Add to the collection
        dfis <- rbind(dfis, dfi)
        
    }
    
    dfis
    
}
