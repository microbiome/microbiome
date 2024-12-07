#' @title Utilities For \code{\link{phyloseq-class}} Slots to Tibbles
#' @description Utility to convert phyloseq slots to tibbles.
#' @details Convert different \code{phyloseq} slots into tibbles.
#' \code{otu_tibble} gets the otu_table in tibble format. 
#' \code{tax_tibble} gets the taxa_table in tibble format. 
#' \code{combine_otu_tax} combines otu_table and taxa_table into one tibble. 
#' @param x \code{\link{phyloseq-class}} object. 
#' @param column.id Provide name for the column which will hold the rownames. 
#'                  of slot.
#' @return A \code{tibble}
#' @examples
#' library(microbiome)
#' data("dietswap")
#' otu_tib <- otu_tibble(dietswap,column.id="FeatureID")
#' tax_tib <- tax_tibble(dietswap,column.id="FeatureID")
#' sample_tib <- sample_tibble(dietswap,column.id="SampleID")
#' otu_tax <- combine_otu_tax(dietswap,column.id = "FeatureID")
#' head(otu_tax)
#' 
#' @name TibbleUtilites
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
NULL

#' @rdname TibbleUtilites
#' @aliases otu_tibble
#' @export
otu_tibble <- function(x, column.id="FeatureID"){
    
    if (any(c("phyloseq", "otu_table") %in% is(x))) {
        
        # Pick OTU matrix
        otu <- abundances(x)
    }
    otu_tibble <- otu %>%  
        as.data.frame() %>% 
        tibble::rownames_to_column(column.id) %>% 
        dplyr::as_tibble()
    return(otu_tibble)
}


#' @rdname TibbleUtilites
#' @aliases tax_tibble
#' @export
tax_tibble <- function(x, column.id="FeatureID"){
    
    if (any(c("phyloseq", "tax_table") %in% is(x))) {
        
        # Pick OTU matrix
        tax <- tax_table(x)
    }
    tax_tibble <- tax %>%  
        as.matrix() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(column.id) %>% 
        dplyr::as_tibble()
    return(tax_tibble)
}

#' @rdname TibbleUtilites
#' @aliases sample_tibble
#' @export
sample_tibble <- function(x, column.id="SampleID"){
    
    if (any(c("phyloseq", "sample_data") %in% is(x))) {
        
        # Pick OTU matrix
        smd <- meta(x)
    }
    sample_tibble <- smd %>% 
        tibble::rownames_to_column(column.id) %>% 
        dplyr::as_tibble()
    return(sample_tibble)
}

#' @rdname TibbleUtilites
#' @aliases combine_otu_tax
#' @export
combine_otu_tax <- function(x, column.id = "FeatureID"){
    otu_tb <- tax_tb <- NULL
    otu_tb <- otu_tibble(x, column.id) 
    tax_tb <- tax_tibble(x, column.id)
    otu_tb <- tax_tb %>% 
        dplyr::left_join(otu_tb, by=column.id)
    return(otu_tb)
}
