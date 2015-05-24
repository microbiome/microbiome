#' get probeset data matrix
#' 
#' 
#'   @param name name
#'   @param level taxonomic level
#'   @param phylogeny.info phylogeny.info
#'   @param probedata oligos vs. samples preprocessed data matrix; 
#'                    absolute scale
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE
#'
#' 
#'   @return probeset data matrix
#'
#' @export
#' @examples 
#'   phylogeny.info <- GetPhylogeny('HITChip', 'filtered')
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   probedata <- read.profiling("frpa", data.dir = data.dir)$oligo
#'   ps <- get.probeset('Akkermansia', 'L2', phylogeny.info, probedata)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get.probeset <- function(name, level, phylogeny.info, probedata, log10 = TRUE) {
    
    phylogeny.info <- as.data.frame(phylogeny.info)

    # Pick probes for this entity
    probes <- retrieve.probesets(phylogeny.info, level, name)
    
    sets <- vector(length = length(probes), mode = "list")
    names(sets) <- names(probes)
    
    for (nam in names(probes)) {
        
        # Pick expression for particular probes (absolute scale)
        p <- intersect(probes[[nam]], rownames(probedata))
        dat <- NULL
        if (length(p) > 0) {
            dat <- probedata[p, ]
            
            dat <- matrix(dat, nrow = length(probes[[nam]]))
            rownames(dat) <- probes[[nam]]
            colnames(dat) <- colnames(probedata)
            
            # Logarithmize probeset?
            if (log10) {
                dat <- log10(dat)
            }
        }
        sets[[nam]] <- dat
        
    }
    
    if (length(sets) == 1) {
        sets <- sets[[1]]
    }
    
    # Return
    sets
    
}