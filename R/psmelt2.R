#' @title Convert \code{\link{phyloseq-class}} object to long data format
#' @description An alternative to psmelt function from 
#'              \code{\link{phyloseq-class}} object.
#' @param x \code{\link{phyloseq-class}} object
#' @param x \code{\link{phyloseq-class}} object
#' @param sample.column A single character string specifying name
#'                      of the column to hold sample names.
#' @param feature.column A single character string specifying name
#'                      of the column to hold OTU or ASV names.
#'                      
#' @examples
#' data("dietswap")
#' ps.melt <- psmelt2(dietswap, sample.column="SampleID", 
#'                    feature.column="Feature") 
#' head(ps.melt)                                         
#' @return A \code{tibble} in long format
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
#' @export                     
psmelt2 <- function(x, sample.column=NULL, feature.column=NULL){
    
    if (class(x)!="phyloseq"){
        stop("Input is not an object of phyloseq class")
    }
    
    if(is.null(sample.column) && any(sample_names(x) %in% sample.column)){
        
        sample.column = "SampleID"
        
    }
    
    if(is.null(sample.column) && !any(sample_names(x) %in% sample.column)){
        
        sample.column = ".SampleID"
        
    }
    
    if(is.null(feature.column)){
        
        feature.column = "FeatureID"
        
    }
    
    otu_tib <- otu_tibble(x, column.id = feature.column)
    tax_tib <- tax_tibble(x, column.id = feature.column)
    sam_tib <- sample_tibble(x, column.id = sample.column)
    otu_tib.m <- otu_tib %>% 
        tidyr::pivot_longer(cols = sample_names(x), names_to = sample.column) %>% 
        dplyr::left_join(tax_tib, by=feature.column) %>% 
        dplyr::left_join(sam_tib, by=sample.column)
    return(otu_tib.m)
}