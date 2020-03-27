#' @title Map Taxonomic Levels
#' @description Map taxa between hierarchy levels.
#' @param taxa taxa to convert; if NULL then considering all taxa in the
#' tax.table
#' @param from convert from taxonomic level 
#' @param to convert to taxonomic level
#' @param data Either a \code{\link{phyloseq}} object or its
#' \code{\link{taxonomyTable-class}} , see the \pkg{phyloseq} package.
#' @return mappings
#' @examples
#' data(dietswap)
#' m <- map_levels('Akkermansia', from='Genus', to='Phylum',
#' tax_table(dietswap))
#' m <- map_levels('Verrucomicrobia', from='Phylum', to='Genus',
#' tax_table(dietswap))
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
map_levels <- function(taxa=NULL, from, to, data) {
    
    if (is.phyloseq(data)) {
        data <- tax_table(data)
    }
    
    # If taxonomy table is from phyloseq, pick the data matrix separately
    if (length(is(data)) == 1 && is(data) == "taxonomyTable") {
        data <- tax_table(data)
    }
    
    df <- data

    if (from == to) {
        df <- list()
        df[[to]] <- factor(taxa)
        df <- as.data.frame(df)
        return(df)
    }
    
    if (is.null(taxa)) {
        taxa <- na.omit(as.character(unique(df[, from])))
    }
    
    # From higher to lower level
    if (length(unique(df[, from])) <= length(unique(df[, to]))) {

        sl <- list()
        for (pt in taxa) {
            
            inds <- which(as.vector(as.character(df[, from])) == pt)
            pi <- df[inds, to]
            sl[[pt]] <- as.character(unique(na.omit(pi)))
            
        }
        
    } else {
        
        # From lower to higher level
        inds <- match(as.character(taxa), df[, from])
        omap <- df[inds, ]
        sl <- omap[, to]
        
    }
    
    as.vector(sl@.Data)
    
}

