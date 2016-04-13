#' @title Individual variability
#' @description Average coefficient of variation (CoV) within individuals for each feature in the data matrix.
#' @param dat Data matrix (samples x features)
#' @param meta Metadata (samples x variables) with the following
#' 	       fields: subject, time
#' @param method Method to calculate the variability within
#'   	  	 individuals: sd (ignores time) or timevar (normalizes by time)
#' @return Vector. Average variability within individuals for each feature.
#' @export
#' @examples
#'   #library(microbiome)
#'   #data(peerj32)
#'   #cv <- variability_individual(peerj32$microbes, peerj32$meta)
#'   #print(head(cv))
variability_individual <- function (dat, meta, method = "timevar") {

  coms <- intersect(rownames(dat), rownames(meta))
  if (length(coms) < 2) {
    stop("Check that the dat and meta arguments have more than 1 samples in common")
  } else {	    
    meta <- meta[coms,]	    
    dat <- dat[coms, ]	    
    subj <- as.character(meta$subject)
    spl <- split(coms, subj) # Sample list for each subject
    spl <- spl[sapply(spl,length) > 1] # only keep subjects with >1 samples

    if (length(spl) == 0) {
      warning("Every subject has only a single sample. Not possible to assess variability within individuals."); 
      return(NA)
    }    

    # For each feature (e.g. microbe) calculate
    # CoV within each subject then average over subjects; time is ignored
    # (intervals are the same for all subjects anyway)
      timevar <- c()
      for (tax in colnames(dat)) {
        stab <- subject_tables(dat[, tax], meta[, c("subject", "time")])

	if (method == "std") {
          m <- sapply(stab, function (tab) {x <- tab$signal; sd(x)/mean(x)})
    	} else if (method == "timevar") {
          m <- sapply(stab, function (tab) {x <- tab$shift[-1]; T <- tab$time[-1]; sqrt(sum((x/T)^2))})
	}
        timevar[[tax]] <- mean(na.omit(m))
      }
  }

  timevar

}
