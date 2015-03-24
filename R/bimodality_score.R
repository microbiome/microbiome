#' Calculate a score based on a given statistical method for detecting bimodality.  
#' The score can be computed using the Hartigan's Dip Test Statistic, coefficient of bimodality, 
#' Silverman's Test or potential analysis. A low score (min = 0) infers unimodality, whereas a
#' high score (max = 1) infers multimodality (or bimodality when using coefficient of bimodality). 
#'
#' @param x Data vector or matrix (samples x phylotypes) 
#' @param method Hartigan's Dip Test Statistic (diptest), coefficient of bimodality (cob),
#' Silverman's Test (silverman) or potential analysis (potential)
#'
#' @return Bimodality score
#'
#' @examples 
#' x <- rnorm(100)
#' bimodality_score(x, method="cob")
#' @export
#' @import devtools 
#' @import earlywarnings
#' @import diptest
#' @import silvermantest
#' @details Hartigan's dip test statistic as described in Hartigan PM
#'          App. Stat, 320-325, 1985.
#'                    
#'          Coefficient of bimodality used in Shade et al. mBio 5(4):e01371-14. 
#'          and picked from Ellison AM Am. J. Bot 74:1280-8, 1987 
#'          
#'          Silverman's bootstrap test for multimodality adopted from
#'          Silverman BW J. Ro. St. Soc 97-99, 1981.
#'           
#'          Potential analysis from R Early Warning Signals Toolbox 
#'          based on Dakos et al. (2012)PL. ONE 7(7):e41010, 2012.
#'          
#' @references 
#'   PM Hartigan, App. Stat, 320-325, 1985.
#'   Shade et al. mBio 5(4):e01371-14, 2014.
#'   AM Ellison, Am. J. Bot 74:1280-8, 1987.
#'   BW Silverman, Jo. R. St. Soc 97-99, 1981. 
#'   Dakos et al. PL. ONE 7(7):e41010, 2012.
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Tineka Blake \email{tineka.blake@@aalto.fi}
#' @keywords utilities
bimodality_score <- function(x, method = "diptest"){
  if (is.vector(x) == TRUE){  
    
    if (method == "diptest") {
      dip <-dip.test(x, simulate.p.value = FALSE, B = 100)
      dip1 <- dip$p.value
      return (1-dip1)
      
    } else if (method == "cob"){
      cob <- coefficient.of.bimodality(x)
      return(cob)
      
    } else if (method == "silverman"){
      sil <- silverman.test(x, 1)
      sil1 <- sil@p_value
      return (1-sil1)
      
    } else if (method == "potential"){
      pa <- multimodality_score(t(x), detection.threshold = 0.9, bs.iterations = 100, detection.limit = 3)
      pa1 <- pa$results$`1`$unimodality.support
      return(1-pa1)
    }
  }
  
  else if ((is.data.frame(x)== TRUE) | (is.matrix(x) == TRUE)){
    if (method == "diptest"){
      dip <- apply(x, 2, function(x) dip.test(x, simulate.p.value = FALSE, B = 100)$statistic)
      return(1-dip)
      
    } else if (method == "cob"){
      cob <- apply(x, 2, coefficient.of.bimodality)
      return (cob)
      
    } else if (method == "silverman"){
      sil <- apply(x, 2, function(x) silverman.test(x,1)@p_value)
      return( 1-sil)
      
    } else if (method == "potential"){
      m <- multimodality_score(t(x), detection.threshold = 0.9, bs.iterations = 100, detection.limit = 3)
      pa <- sapply(m$results, function (x){x$unimodality.support})
      return(1-pa)
    }
  } 
} 