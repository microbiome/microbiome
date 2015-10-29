#' @title intermediate_stability
#' @description Quantify intermediate stability with respect to a given
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


