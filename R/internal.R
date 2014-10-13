
get.file.method <- function( f ) {
    
    method <- NULL
    
    if (f == "rpa") {
        method <- "rpa"
    }
    if (f == "frpa") {
        method <- "frpa"
    }
    if (!length(grep("sum", f)) == 0) {
        method <- "sum"
    }
    if (!length(grep("ave", f)) == 0) {
        method <- "ave"
    }
    if (!length(grep("nmf", f)) == 0) {
        method <- "nmf"
    }
    
    method
    
}

get.file.level <- function(f) {
    
    level <- NULL
    
    if (!length(grep("oligo", f)) == 0) {
        level <- "oligo"
    }
    if (!length(grep("species", f)) == 0) {
        level <- "species"
    }
    if (!length(grep("L0", f)) == 0) {
        level <- "L0"
    }
    if (!length(grep("L1", f)) == 0) {
        level <- "L1"
    }
    if (!length(grep("L2", f)) == 0) {
        level <- "L2"
    }
    if (!length(grep("phylogeny", f)) == 0) {
        level <- "phylogeny.info"
    }
    
    level
    
}

#' Species-levels mappings
#'
#' @param spec species
#' @param phylogeny.info phylogeny.info
#'
#' @return species-levels mappings
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal

species2levels <- function(spec, phylogeny.info) {
    
    phylogeny.info <- polish.phylogeny.info(phylogeny.info)
    
    omap <- phylogeny.info[match(as.character(spec), 
               phylogeny.info$species), c("species", 
        "L2", "L1")]
    omap[["species"]] <- factor(omap[["species"]])
    omap[["L1"]] <- factor(omap[["L1"]])
    omap[["L2"]] <- factor(omap[["L2"]])
    
    omap
}

#' level2-level1 mappings
#'
#' @param l2 level2 phylotypes
#' @param phylogeny.info phylogeny.info
#'
#' @return mappings
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal

level2TOlevel1 <- function(l2, phylogeny.info) {
    
    phylogeny.info <- polish.phylogeny.info(phylogeny.info)
    
    omap <- phylogeny.info[match(as.character(l2), 
               phylogeny.info[["L2"]]), c("L2", 
        "L1")]
    omap[["L2"]] <- factor(omap[["L2"]])
    omap[["L1"]] <- factor(omap[["L1"]])
    
    omap
}


#' L2-L0 mappings
#'
#' @param l2 level2 phylotypes
#' @param phylogeny.info phylogeny.info
#'
#' @return mappings
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal

level2TOlevel0 <- function(l2, phylogeny.info) {
    
    phylogeny.info <- polish.phylogeny.info(phylogeny.info)
    
    omap <- phylogeny.info[match(as.character(l2), 
              phylogeny.info[["L2"]]), c("L2", 
        "L0")]
    omap[["L2"]] <- factor(omap[["L2"]])
    omap[["L0"]] <- factor(omap[["L0"]])
    
    omap
}


#' L1-L0 mappings
#'
#' @param l1 level1 phylotypes
#' @param phylogeny.info phylogeny.info
#'
#' @return mappings
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal

level1TOlevel0 <- function(l1, phylogeny.info) {
    
    phylogeny.info <- polish.phylogeny.info(phylogeny.info)
    
    omap <- phylogeny.info[match(as.character(l1), 
                phylogeny.info[["L1"]]), c("L1", 
        "L0")]
    omap[["L1"]] <- factor(omap[["L1"]])
    omap[["L0"]] <- factor(omap[["L0"]])
    
    omap
}



#' oligo-species mappings
#'
#' @param oligo oligos
#' @param phylogeny.info phylogeny.info
#'
#' @return mappings
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal

oligoTOhigher <- function(oligo, phylogeny.info, level.to = "species") {
    
    phylogeny.info <- polish.phylogeny.info(phylogeny.info)
    
    oligo.species <- list()
    for (oli in as.character(oligo)) {
        oligo.species[[oli]] <- sort(unique(
              phylogeny.info[which(phylogeny.info[["oligoID"]] == 
            oli), level.to]))
    }
    
    oligo.species
} 
