#' @title Import HITChip data
#' @description Import HITChip output into phyloseq format.
#' @param data.dir Profiling script output directory for reading the data. 
#' @param method Probe summarization method ("rpa" or "sum")
#' @param detection.threshold Taxon absence/presence thresholds (typically 10^1.8 for HITChip)
#' @param verbose verbose
#' @return data matrix (phylo x samples)
#' @export
#' @examples
#'  \dontrun{
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   dat <- import_hitchip(data.dir)
#' }
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
import_hitchip <- function(data.dir, method = "frpa", detection.threshold = 0, verbose = F) {

  # Read	     
  if ( verbose ) { message(paste("Reading Chip data from", data.dir)) }

  res <- list()

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  res[["probedata"]] <- tab

  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  if (!file.exists(f)) {
    # Old outputs had this name
    f <- paste(data.dir, "/phylogeny.filtered.tab", sep = "")    
  }
  taxonomy <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  res[["taxonomy"]] <- taxonomy
  
  # Read unfiltered taxonomy table
  f <- paste(data.dir, "/taxonomy.full.tab", sep = "")
  if (!file.exists(f)) {
    # Old outputs had this name
    f <- paste(data.dir, "/phylogeny.full.tab", sep = "")    
  }
  taxonomy.full <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  res[["taxonomy.full"]] <- taxonomy.full

  # Read sample metadata      
  f <- paste(data.dir, "/meta.tab", sep = "")
  if (file.exists(f)) {
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    rownames(tab) <- tab$sample
    meta <- tab
    res[["meta"]] <- meta
  }

  # -----------------------------------

  # Only pick probe-level data, taxonomy and metadata
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


