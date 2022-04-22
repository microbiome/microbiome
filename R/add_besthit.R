#' @title Adds \code{best_hist} to a \code{\link{phyloseq-class}} Object
#' @description Add the lowest classification for an OTU or ASV.
#' @details Most commonly it is observed that taxa names are either OTU ids or 
#'          ASV ids. In such cases it is useful to know the taxonomic identity.
#'          For this purpose, \code{best_hist} identifies the best available 
#'          taxonomic identity and adds it to the OTU ids or ASV ids. If genus 
#'          and species columns are present in input the function internally 
#'          combines the names.  
#' @param x \code{\link{phyloseq-class}} object
#' @param sep separator e.g. ASV161:Roseburia
#' @return  \code{\link{phyloseq-class}} object \code{\link{phyloseq-class}}
#' @export
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiome)
#' data(dietswap)
#' p0.f <- add_besthit(atlas1006, sep=":")
#' }
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities

add_besthit <- function(x, sep=":"){
  
  Class<-Domain<- Family<- Genus<- Genus.Species<- NULL
  Order<- Phylum<- Species<-NULL
  
  x.nw <- x
  if(length(rank_names(x.nw))== 6){
    colnames(tax_table(x.nw)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  }
  if(length(rank_names(x.nw))==7){
    colnames(tax_table(x.nw)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  
  tax.tib <- .get_taxa_tib_unite(x)
  
  tax.tib <- tax.tib %>% 
    dplyr::mutate(Domain =ifelse(is.na(Domain), "Unclassifed", Domain),
                  Phylum =ifelse(is.na(Phylum), Domain, Phylum),
                  Class =ifelse(is.na(Class), Phylum, Class),
                  Order =ifelse(is.na(Order), Class, Order),
                  Family =ifelse(is.na(Family), Order, Family),
                  Genus =ifelse(is.na(Genus), Family, Genus)) 
  if(length(rank_names(x))==7){
    tax.tib <- tax.tib %>%
      dplyr::mutate(Species =ifelse(is.na(Species), Genus, Species))
  }
  
  best_hit <- paste0(taxa_names(x), sep,tax.tib[,ncol(tax.tib)])
  
  taxa_names(x) <- best_hit
  return(x)
}



.get_taxa_tib_unite <- function(x){
  
  Genus<- Species <- Genus.Species<- NULL
  tax.tib <- tax_table(x) %>% 
    as.matrix() %>% 
    as.data.frame() 
  
  #n.rk <- length(rank_names(x))
  if(any(rank_names(x) == "Species") && any(rank_names(x) == "Genus")){
    
    tax.tib <- tax.tib %>% 
      dplyr::mutate(Genus.Species = ifelse(!is.na(Species), 
                                           paste0(Genus, ".", Species), Species)) %>%
      dplyr::select(-Species) %>%
      dplyr::rename(Species = Genus.Species)
    
  }
  return(tax.tib)
}



