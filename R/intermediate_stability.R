#' @title Intermediate stability
#' @description Quantify intermediate stability with respect to a given reference point. 
#' @param x \pkg{phyloseq} object.
#'          Includes otu_table (variables x samples) and
#' 	    sample_data data.frame (samples x features) with 'subject'
#'	    and 'time' field for each sample.
#' @param reference.point Calculate stability of the data w.r.t. this point. By default the intermediate range is used (min + (max - min)/2). If a vector of points is provided, then the scores will be calculated for every point and a data.frame is returned.
#' @param method 'lm' (linear model) or 'correlation';
#'               the linear model takes time into account as a covariate 
#' @param output Specify the return mode. Either the "full" set of stability
#'        analysis outputs, or the "scores" of intermediate stability.
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
#' data, and subjects with less than two time point are excluded. The absolute
#' count data x is logarithmized before the analysis with the log10(1 + x)
#' trick to circumvent logarithmization of zeroes.
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
intermediate_stability <- function (x, reference.point = NULL, method = "correlation", output = "scores") {

  if (length(reference.point) > 1 && output == "scores") {
    scores <- c()
    for (i in 1:length(reference.point)) {
      scores[[i]] <- intermediate_stability(x, reference.point = reference.point[[i]], method = method, output = output)
     }

     if (ntaxa(x) > 1) {
       scores <- as.data.frame(scores, nrow = ntaxa(x))
     } else {
       scores <- as.data.frame(t(as.matrix(scores, ncol = length(reference.point))))
     }

     colnames(scores) <- as.character(reference.point)
     rownames(scores) <- taxa_names(x)

     return(scores)
  }
  

  # Logarithmize the data with log10(1 + x) trick
  pseq <- x		       
  x <- log10(1 + t(otu_table(pseq)@.Data))
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


