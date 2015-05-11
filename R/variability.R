#' Variability analysis. Calculates average Pearson '
#' correlation between samples in the input data and picks the lower '
#' triangular matrix to avoid duplicating the correlations. Returns 
#' correlations and stability estimate (average of the correlations). 
#' Can also be used to calculate stability between two data sets. 
#' Then provide two data sets as inputs.
#'
#' @param dat1 data matrix samples vs phylotypes (in log10 scale)
#' @param dat2 Optional. Second data matrix samples vs. phylotypes. 
#'          Provide this to calculate stability between two (paired) 
#'          data sets.
#' @param method Correlation method (see ?cor)
#'
#' @return List with correlations and astability estimate
#' @import dplyr
#'
#' @export
#' @examples 
#'   library(microbiome)
#'   data.peerj32 <- download_microbiome("peerj32")
#'   x <- data.peerj32$microbes
#'   s <- estimate_variability(x[, 1:5])
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
estimate_variability <- function(x, meta, type, group_by = "group", method = "spearman") {
    
    # Split the data by group
    if (!group_by %in% names(meta)) {
      meta[[group_by]] <- rep("completedata", nrow(meta))
    }
    datasets <- split(x, droplevels(meta[[group_by]]))

    if (type == "interindividual") {
      
      # Within-matrix stability NOTE: earlier this was calculated as
      # the average of upper triangular correlation matrix This is
      # heavily biased since the values are dependent Now replaced
      # by calculating correlations against the mean of the whole
      # sample set cors <- lower.triangle(cor(dat1))
      dfs <- NULL      
      for (ds in names(datasets)) {
        dat1 <- datasets[[ds]]
        cors <- as.vector(cor(t(dat1), matrix(colMeans(dat1)), method = method))
        dfs <- rbind(dfs, data.frame(group = rep(ds, length(cors)),
                   	             sample = rownames(dat1),
  		   		     correlation = cors))
	#variability[[ds]] <- list(correlations = cors, variability = mean(cors))
      }
  
      pval <- anova(lm(correlation ~ group, data = dfs))[["Pr(>F)"]][[1]]
      stats <- dfs %>% group_by(group) %>% summarize(mean = mean(correlation), sd = sd(correlation))
      variability <- list(data = dfs, statistics = stats, p.value = pval)

    } else if (type == "intraindividual") {

      variability <- list()
      dfs <- NULL
      for (ds in names(datasets)) {
      
        # Pick the data and metadata for this group
        xsub <- datasets[[ds]]        
	msub <- meta[rownames(xsub),]

        # Use interindividual functionality to assess correlations
	# within subjects. Subjects are used as groups
        datasets2 <- split(xsub, droplevels(msub[["subject"]]))
	cors <- c()
        for (subj in names(datasets2)) {
          dats <- datasets2[[subj]]
          cors[[subj]] <- cor(unlist(dats[1,]), unlist(dats[2,]), method = method)
        }

        dfs <- rbind(dfs, data.frame(group = rep(ds, length(cors)),
                   	             subject = names(cors),
  		   		     correlation = cors))

      }

      # Between time point correlations within subjects
      # and the mean over those correlations
      pval <- anova(lm(correlation ~ group, data = dfs))[["Pr(>F)"]][[1]]
      stats <- dfs %>% group_by(group) %>% summarize(mean = mean(correlation), sd = sd(correlation))
      variability <- list(data = dfs, statistics = stats, p.value = pval)

    }
    
    variability
    
}



#' variability_individual
#' 
#' Average CoV within individuals for each feature in the data matrix
#'
#' @param dat Data matrix (samples x features)
#' @param meta Metadata (samples x variables) with the following
#' 	       fields: subject, time
#' @param method Method to calculate the variability within
#'   	  	 individuals: sd (ignores time) or timevar (normalizes by time)
#'
#' 
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


subject_tables <- function (x, meta) {

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




