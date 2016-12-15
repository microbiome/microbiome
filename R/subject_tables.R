subject_tables <- function (x, meta) {

  stop("TODO: replaced now with timesort_subjects.R")

  # Focus on the signal from specific taxon
  meta$signal <- x

  # Pick data for each subject separately
  spl <- split(meta, as.character(meta$subject))

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




