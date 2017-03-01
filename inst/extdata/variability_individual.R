#' @title Individual Variability
#' @description Average coefficient of variation (CoV) within individuals
#'              for each taxa.
#' @param x phyloseq object
#' @param method "CoV" (coef. of variation ie. std/mean) or
#'               "std" (standard deviation) or
#'               "timevar" (average shift normalized by time)
#' @details The sample metadata should contain the subject and time fields
#' @return A vector of average variability within individuals for each taxa.
#' @export
#' @examples
#'   library(microbiome)
#'   data(atlas1006)
#'   cv <- variability_individual(atlas1006)
variability_individual <- function (x, method = "CoV") {

  # TODO speed up calculations

  dat <- t(abundances(x))
  meta <- sample_data(x)		         

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

    # For each taxa calculate
    timevar <- c()
    for (tax in colnames(dat)) {
      m <- meta[, c("subject", "time")]
      m$signal <- dat[, tax]
      stab <- time_sort(m)

      if (method == "CoV") {	
          # CoV within each subject then average over subjects; time is ignored
          m <- sapply(stab, function (tab) {x <- tab$signal; sd(x)/mean(x)})
      } else if (method == "std") {	
          # CoV within each subject then average over subjects; time is ignored
          m <- sapply(stab, function (tab) {sd(tab$signal)})
      } else if (method == "timevar") {
	  # Change normalized by time
          m <- sapply(stab, function (tab) {x <- tab$shift[-1]; T <- tab$time[-1]; sqrt(sum((x/T)^2))})
      }
      
      timevar[[tax]] <- mean(na.omit(m))
    }

  }

  timevar

}
