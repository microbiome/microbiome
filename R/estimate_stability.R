#' Quantify intermediate stability with respect to a given
#' reference point. 
#'
#' @param x \pkg{phyloseq} object.
#'          Includes otu_table (variables x samples) and
#' 	    sample_data data.frame (samples x features) with 'subject'
#'	    and 'time' field for each sample.
#'           
#' @param reference.point Optional. Calculate stability of the
#'                        data w.r.t. this point. By default the
#'                        intermediate range is used (min + (max - min)/2)
#' @param method 'lm' (linear model) or 'correlation';
#'               the linear model takes time into account as a covariate 
#' @param output Specify the return mode. Either the "full" set of stability
#'        analysis outputs, or the "scores" of intermediate stability.
#' 
#' @return A list with following elements: 
#' 	     stability: estimated stability
#'	     data: processed data set used in calculations	    
#'
#' @details Decomposes each column in x into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' "correlation" method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The "lm" method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance). Samples with missing
#' data, and subjects with less than two time point are excluded.	   
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples
#' # Example data
#' #library(microbiome)
#' #x <- download_microbiome("atlas1006")
#  #x <- prune_taxa(c("Akkermansia", "Dialister"), x)
#' #res <- intermediate_stability(x, reference.point = NULL)
#' #s <- sapply(res, function (x) {x$stability})
#' @keywords utilities
intermediate_stability <- function (x, reference.point = NULL, method = "correlation", output = "full") {

  # Logarithmize the data
  pseq <- x		       
  x <- log10(t(otu_table(pseq)@.Data))
  meta <- sample_data(pseq)

  # Estimate stabilities for each OTU
  stability <- list()
  for (tax in colnames(x)) {
    df <- meta
    df$data <- x[, tax]
    stability[[tax]] <- estimate_stability(df, 
    		     	  reference.point = reference.point, 
		     	  method = method)
  }

  if (output == "full") {
    return(stability)
  } else if (output == "scores") {
    intermediate.stability <- sapply(stability, function (x) {x$stability})
    return(intermediate.stability)  
  }

}


#' estimate_stability
#'
#' Description: Quantify intermediate stability with respect to a given reference point. 
#'
#' @param df Combined input data vector (samples x variables) and metadata data.frame (samples x features)
#'           with the 'data', 'subject' and 'time' field for each sample 
#'           
#' @param reference.point Optional. Calculate stability of the data w.r.t. this point. 
#'                        By default the intermediate range is used (min + (max - min)/2)
#' @param method "lm" (linear model) or "correlation"; the linear model takes time into account 
#' 	         as a covariate 
#' 
#' @return A list with following elements: 
#' 	     stability: estimated stability
#'	     data: processed data set used in calculations	    
#'
#' @details Decomposes each column in x into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' "correlation" method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The "lm" method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance). 
#' Samples with missing data, and subjects with less than two time point are excluded.	   
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples 
#'   #df <- data.frame(list(
#'   #	  subject = rep(paste("subject", 1:50, sep = "-"), each = 2), 
#'   #	  time = rep(1:2, 50),
#'   #	  data = rnorm(100)))
#'   # s <- estimate_stability_single(df, reference.point = NULL, method = "lm")
#'
#' @keywords utilities

estimate_stability <- function (df, reference.point = NULL, method = "lm") {

  # Remove NAs
  df <- df[!is.na(df$data),]

  # Detect intermediate value in the overall data if reference point not given
  if (is.null(reference.point)) {
    reference.point <- mean(range(df$data))
  }

  # Remove subjects with only one measurement
  df <- df[df$subject %in% names(which(table(df$subject) > 1)),]

  if (nrow(df) < 2) {warning("No subjects with time series in estimate_stability. Returninng NULL"); return(NULL)} 

  # Split data by subject
  spl <- split(df, as.character(df$subject))

  dfis <- NULL

  for (spli in spl) {

    # Ensure the measurements are ordered in time
    spli <- spli[order(spli$time),]

    # Calculate differences between consecutive time points; 
    # and the initial values; and their distance from reference
    data.difs <- diff(spli$data)
    time.difs <- diff(spli$time)
    start.points <- spli$data[-nrow(spli)]

    start.reference.distance <- start.points - reference.point

    # Organize into data frame
    dfi <- data.frame(change = data.difs, time = time.difs, start = start.points, start.reference.distance = start.reference.distance)

    # Add to the collection
    dfis <- rbind(dfis, dfi)

  }

  dfis.left <- subset(dfis, start.reference.distance < 0)
  dfis.right <- subset(dfis, start.reference.distance > 0)

  # Simplified stability calculation (do not consider time effect)
  stability <- NULL
  stability.left <- stability.right <- NA
  if (method == "correlation") {
    # For each subject, check distance from the stability point
    # at the baseline time point
    baseline.distance <- abs(dfis$start.reference.distance)
    # For each subject, calculate deviation between the first and second time point
    followup.distance <- abs(dfis$change)
    stability <- cor(baseline.distance, followup.distance)  

    if (nrow(dfis.left)>10){ 
      # Negative values for low stability
      baseline.distance <- abs(dfis.left$start.reference.distance)
      followup.distance <- dfis.left$change
      stability.left <- cor(baseline.distance, followup.distance)  

    }

    if (nrow(dfis.right)>10) {
      baseline.distance <- abs(dfis.right$start.reference.distance)
      followup.distance <- dfis.right$change
      stability.right <- -cor(baseline.distance, followup.distance)  
    }


  } else if (method == "lm") {
    # Advanced calculation, take time into account with linear model (also possible to check p-values later if needed)
    stability <- coef(summary(lm(abs(change) ~ time + abs(start.reference.distance), data = dfis)))["abs(start.reference.distance)", "Estimate"]
  
    if (nrow(dfis.left)>10){ 
      # Negative values for low stability
      stability.left <- coef(summary(lm(change ~ time + abs(start.reference.distance), data = dfis.left)))["abs(start.reference.distance)", "Estimate"]
    }
    if (nrow(dfis.right)>10){ 
      # Negative values for low stability
      stability.right <- -coef(summary(lm(change ~ time + abs(start.reference.distance), data = dfis.right)))["abs(start.reference.distance)", "Estimate"]

    }
  }

  list(stability = stability, stability.right = stability.right, stability.left = stability.left, data = dfis)

}
