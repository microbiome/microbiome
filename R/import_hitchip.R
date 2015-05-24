#' Import HITChip output into phyloseq format
#'
#' @param data.dir Profiling script output directory for reading the data. 
#' @param method Probe summarization method ("rpa" or "sum")
#' @param detection.threshold Taxon absence/presence thresholds (typically 10^1.8 for HITChip)
#' 
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples 
#'   # data.dir <- system.file("extdata", package = "microbiome")
#'   # dat <- import_hitchip(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
import_hitchip <- function(data.dir, method = "frpa", detection.threshold = 0) {

  # Read	     
  res <- read.profiling(data.dir, method = method)
  taxonomy <- res$taxonomy
  probedata <- res$probedata
  meta <- res$meta

  # Summarize probes into abundance table
  level <- "species"
  abu <- summarize_probedata(probedata = probedata,
      	 	             taxonomy = taxonomy,
      	 		     level = level,
			     method = method)

  # Convert the object into phyloseq format
  levels <- intersect(c("L0", "L1", "L2", "species"), colnames(taxonomy))
  taxonomy <- unique(taxonomy[, levels])
  rownames(taxonomy) <- taxonomy[, "species"]
  coms <- intersect(rownames(taxonomy), rownames(abu))
  abu <- abu[coms,]
  taxonomy <- taxonomy[coms,]

  pseq <- hitchip2physeq(t(abu), meta, taxonomy, detection.limit = detection.threshold)

  return(pseq) 
    
} 

