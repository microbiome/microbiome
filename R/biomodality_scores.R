#' Calculate Biomodality scores 
#'
#' @details Calculate bimodality scores using, Hartigan's Dip Test Statistic, 
#' Coeffefient of bimodality, Silvermans Test and Potential analysis.
#'
#' @param x Data (vector or matrix)
#' @param method Statistical test 
#'
#' @return bimodality score
#'
#' @examples 
#'   library(microbiome)
#'   bimodality_score(x, method= "bc")
#'
#' @export
#' @references 
#'   Code modified from the original source:
#'   r-bloggers.com/measuring-associations-between-non-numeric-variables/
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords bimodality

bimodality_score <- function(x, method = ...){
  if (is.vector(x) == TRUE){  
    
    if (method == "diptest") {
      dip.test(x, simulate.p.value = FALSE, B = 100)
      
    } else if (method == "bc"){
      cof <- coefficient.of.bimodality(x)
      return (1-cof)
      
    } else if (method == "silverman"){
      silverman.test(x, 1)
      
    } else if (method == "potential"){
      multimodality_score(t(x), detection.threshold = 0.9, bs.iterations = 100, detection.limit = 3)
    }
  }
  
  else if ((is.data.frame(x)== TRUE) | (is.matrix(x) == TRUE)){
    if (method == "diptest") {
      apply(x, 2, function (x) {bimodality(x, method = "diptest")$p.value})
      
    } else if (method == "bc"){
      bs <- apply(x, 2, coefficient.of.bimodality)
      return (1-bs)
      
    } else if (method == "silverman"){
      apply(x, 2, function (x) {bimodality(x, method = "silverman")@p_value})
      
    } else if (method == "potential"){
      m <- multimodality_score(t(x), detection.threshold = 0.9, bs.iterations = 100, detection.limit = 3)
      sa <- sapply(m$results, function (x){x$unimodality.support})
    }
  } 
}
