#' get probeset data matrix
#' 
#' 
#'   @param name name
#'   @param level taxonomic level
#'   @param taxonomy taxonomy
#'   @param probedata oligos vs. samples preprocessed data matrix; 
#'                    absolute scale
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE
#'
#' 
#'   @return probeset data matrix
#'
#' @export
#' @examples 
#'   taxonomy <- GetPhylogeny('HITChip', 'filtered')
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   probedata <- read.profiling(data.dir, "rpa")$probedata
#'   ps <- get.probeset('Akkermansia', 'L2', taxonomy, probedata)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get.probeset <- function(name, level, taxonomy, probedata, log10 = TRUE) {
    
    taxonomy <- as.data.frame(taxonomy)

    # Pick probes for this entity
    probes <- retrieve.probesets(taxonomy, level, name)
    
    sets <- vector(length = length(probes), mode = "list")
    names(sets) <- names(probes)
    
    for (nam in names(probes)) {
        
        # Pick expression for particular probes (absolute scale)
        p <- intersect(probes[[nam]], rownames(probedata))
        dat <- NULL
        if (length(p) > 0) {
            dat <- probedata[p, ]  
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