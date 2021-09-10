#' @title Richness Index
#' @description Community richness index.
#' @inheritParams alpha
#' @param index "observed" or "chao1"
#' @param detection Detection threshold. Used for the "observed" index.
#' @return A vector of richness indices
#' @details By default, returns the richness for multiple detection thresholds
#' defined by the data quantiles. If the detection argument is provided,
#' returns richness with that detection threshold. The "observed" richness
#' corresponds to index="observed", detection=0.
#' @export
#' @examples
#' data(dietswap)
#' d <- richness(dietswap, detection=0)
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso alpha
#' @keywords utilities
richness <- function(x, index = c("observed", "chao1"), detection=0) {

    index <- gsub("richness_", "", index)

    # This already ensures that taxa are on the rows     
    x <- abundances(x)
    tab <- NULL
    index <- tolower(index)

    if ("observed" %in% index) {
        tab <- richness_help(x, detection)
        if (is.vector(tab)) {
            tab <- as.matrix(tab, ncol=1)
            colnames(tab) <- gsub("%", "", as.character(detection))
        }
    }

    if ("chao1" %in% index) {
        chao <- chao1(x)
        tab <- cbind(tab, chao1 = chao)
    }

    colnames(tab) <- gsub("^richness_0$", "observed", colnames(tab))
    colnames(tab) <- gsub("^0$", "observed", colnames(tab))    

    as.data.frame(tab)

    
}


richness_help <- function(x, detection=NULL) {
    
    # Pick data
    otu <- x
    
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


# This calculates Chao1 per sample
chao1 <- function (x) {

    x <- as.matrix(x)

    # This can be done to calculate Chao1 for the total sample
    #x <- rowSums(x)
    apply(x, 2, function (xi) {chao1_per_sample(xi)})
    
}


chao1_per_sample <- function (x) {

    s0 <- sum(x>0, na.rm = TRUE)
    s1 <- sum(x==1, na.rm = TRUE)
    s2 <- sum(x==2, na.rm = TRUE)
    # if ((s1-s2)^2==(s1+s2)^2) {
        # Bias corrected version / See issue #150
        r <- s0+s1*(s1-1)/((s2+1)*2)
    #} else {
    #    r <- s0+s1^2/(s2*2)
    #}

    r

}

