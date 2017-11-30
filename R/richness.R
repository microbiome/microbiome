#' @title Richness Index
#' @description Community richness index.
#' @inheritParams global
#' @param detection Detection threshold.
#' @return A vector of richness indices
#' @details By default, returns the richness for multiple detection thresholds
#' defined by the data quantiles. If the detection argument is provided,
#' returns richness with that detection threshold.
#' @export
#' @examples
#' data(dietswap)
#' d <- richness(dietswap, detection=0)
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso global
#' @keywords utilities
richness <- function(x, detection = NULL) {
        
    tab <- richness_help(x, detection)

    if (is.vector(tab)) {
        tab <- as.matrix(tab, ncol=1)
        colnames(tab) <- gsub("%", "", as.character(detection))
    }

    as.data.frame(tab)
    
}


richness_help <- function(x, detection=NULL) {
    
    # Pick data
    otu <- abundances(x)
    
    # Check with varying detection thresholds
    if (is.null(detection) || length(detection) > 0) {
        
        if (is.null(detection)) {
            
            ths <- quantile(as.vector(otu), c(0, 0.2, 0.5, 0.8))
            
        } else {
            
            ths <- detection
            names(ths) <- as.character(ths)
            
        }
        
        tab <- NULL
        for (th in ths) {
            r <- colSums(otu > th)
            tab <- cbind(tab, r)
        }
        
        colnames(tab) <- gsub("%", "", names(ths))
        r <- tab
        
    } else {
        
        r <- colSums(otu > detection)
        names(r) <- colnames(otu)
        
    }
    
    r
    
}


