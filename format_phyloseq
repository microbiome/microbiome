#' @title Formatting the Phyloseq Object
#' @description Format the phyloseq object to correct the taxonomy in phyloseq object (tax_table).
#' @details Most commonly it is observed that the taxonomy file has classification until a given taxonomic level.
#'          Hence, to avoid loss of OTU information while using the function tax_glom() for merging at a specific taxonomic level.
#'          we will fill the empty cells with the maximum classification available along with the OTU number. This code is a slight modification.
#'          the code from  \pkg{ampvis} \code{\link{phyloseq-class}}. Here, we directly take the phyloseq object as input and make the necessary formatting.
#' @param x \code{\link{phyloseq-class}} object
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data(DynamicsIBD)
#'     p0 <- DynamicsIBD
#'     p0.f <- format_phyloseq(p0)
#'           }
#' @keywords utilities
format_phyloseq <- function(x)
{
  if (ncol(tax_table(x)) == 6)
  {
    pobj <- x
    colnames(tax_table(pobj)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
    x <- as.data.frame(pobj@tax_table)
    x$Domain <- gsub("k__", "", x$Domain)
    x$Phylum <- gsub("p__", "", x$Phylum)
    x$Phylum <- gsub("c__", "", x$Phylum)
    x$Class <- gsub("c__", "", x$Class)
    x$Order <- gsub("o__", "", x$Order)
    x$Family <- gsub("f__", "", x$Family)
    x$Genus <- gsub("g__", "", x$Genus)
    tax <- mutate(x, Domain, Domain = ifelse(Domain == "", "Unclassified", Domain)) %>% mutate(Phylum,
                                                                                               Phylum = ifelse(Phylum == "", paste("k__", Domain, "_", rownames(x), sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__",
                                                                                          Phylum, "_", rownames(x), sep = "")), Class)) %>% mutate(Order, Order = ifelse(Order ==
                                                                                                                                                                           "", ifelse(grepl("__", Class), Class, paste("c__", Class, "_", rownames(x), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__",
                                                                                           Order, "_", rownames(x), sep = "")), Family)) %>% mutate(Genus, Genus = ifelse(Genus ==
                                                                                                                                                                            "", ifelse(grepl("__", Family), Family, paste("f__", Family, "_", rownames(x), sep = "")),
                                                                                                                                                                          Genus))
    me <- as.data.frame(pobj@tax_table)
    me$Domain <- tax$Domain
    me$Phylum <- tax$Phylum
    me$Class <- tax$Class
    me$Order <- tax$Order
    me$Family <- tax$Family
    me$Genus <- tax$Genus
    taxmat <- as.matrix(me)
    new.tax <- tax_table(taxmat)
    tax_table(pobj) <- new.tax
    return(pobj)
  } else if (ncol(tax_table(x)) == 7)
  {
    pobj2 <- x
    colnames(tax_table(pobj2)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    x <- as.data.frame(pobj2@tax_table)
    x$Domain <- gsub("k__", "", x$Domain)
    x$Phylum <- gsub("p__", "", x$Phylum)
    x$Phylum <- gsub("c__", "", x$Phylum)
    x$Class <- gsub("c__", "", x$Class)
    x$Order <- gsub("o__", "", x$Order)
    x$Family <- gsub("f__", "", x$Family)
    x$Genus <- gsub("g__", "", x$Genus)
    x$Species <- gsub("s__", "", x$Species)
    tax <- mutate(x, Domain, Domain = ifelse(Domain == "", "Unclassified", Domain)) %>% mutate(Phylum,
                                                                                               Phylum = ifelse(Phylum == "", paste("k__", Domain, "_", rownames(x), sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__",
                                                                                          Phylum, "_", rownames(x), sep = "")), Class)) %>% mutate(Order, Order = ifelse(Order ==
                                                                                                                                                                           "", ifelse(grepl("__", Class), Class, paste("c__", Class, "_", rownames(x), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__",
                                                                                           Order, "_", rownames(x), sep = "")), Family)) %>% mutate(Genus, Genus = ifelse(Genus ==
                                                                                                                                                                            "", ifelse(grepl("__", Family), Family, paste("f__", Family, "_", rownames(x), sep = "")),
                                                                                                                                                                          Genus)) %>% mutate(Species, Species = ifelse(Species == "", ifelse(grepl("__", Genus), Genus,
                                                                                                                                                                                                                                             paste("f__", Genus, "_", rownames(x), sep = "")), Species))
    me <- as.data.frame(pobj2@tax_table)
    me$Domain <- tax$Domain
    me$Phylum <- tax$Phylum
    me$Class <- tax$Class
    me$Order <- tax$Order
    me$Family <- tax$Family
    me$Genus <- tax$Genus
    me$Species <- tax$Species
    taxmat <- as.matrix(me)
    new.tax <- tax_table(taxmat)
    tax_table(pobj2) <- new.tax
  }
  return(pobj2)
}
