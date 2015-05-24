#' Read HITChip run.profiling.script output into R
#'
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
#' @param method Select the preprocessing method that you like to check
#' 
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples 
#'   # data.dir <- system.file("extdata", package = "microbiome")
#'   # dat <- read.profiling(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read.profiling <- function(data.dir, method, verbose = TRUE) {

  if (verbose) {message(paste("Reading Chip data from", data.dir))}

  results <- list()

# Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  results[["probedata"]] <- tab

  # Read abundance tables
  for (s in c("L0", "L1", "L2", "species")) {
    f <- paste(data.dir, "/", s, "-", method, ".tab", sep = "")
    if (file.exists(f)) {
      tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
      colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
      results[[s]] <- tab
    }
  }
  
  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  if (!file.exists(f)) {
    # Old outputs had this name
    f <- paste(data.dir, "/phylogeny.filtered.tab", sep = "")    
  }
  taxonomy <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  results[["taxonomy"]] <- taxonomy
  
  # Read unfiltered taxonomy table
  f <- paste(data.dir, "/taxonomy.full.tab", sep = "")
  if (!file.exists(f)) {
    # Old outputs had this name
    f <- paste(data.dir, "/phylogeny.full.tab", sep = "")    
  }
  taxonomy.full <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  results[["taxonomy.full"]] <- taxonomy.full

  # Read sample metadata      
  f <- paste(data.dir, "/meta.tab", sep = "")
  if (file.exists(f)) {
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    rownames(tab) <- tab$sample
    meta <- tab
    results[["meta"]] <- meta
  }

  results
    
} 


#' Read HITChip output and preprocess into phyloseq format
#'
#' @param data.dir Profiling script output directory for reading the data. 
#' @param method Probe summarization method ("rpa", "frpa", or "sum")
#' @param detection.threshold Taxon absence/presence thresholds (typically 10^1.8 for HITChip)
#' 
#' Converts the probe-level data matrix and probe-level taxonomy table to phylotype-level 
#' HITChip data. Returns the probe-level data (data matrix and taxonomies) and the phylotype-level 
#' phyloseq object. There are two versions of probe-level taxonomy. The full version includes all probes 
#' in the probe-level data. The filtered version includes those probes that have been used to aggregate
#' probes into phylotype level.
#'
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples 
#'   # data.dir <- system.file("extdata", package = "microbiome")
#'   # dat <- read_hitchip(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read_hitchip <- function(data.dir, method = "rpa", detection.threshold = 0) {

  # Convert to phyloseq format. Includes summarization to species level.
  pseq <- import_hitchip(data.dir, method = method, detection.threshold = detection.threshold)

  # Read data
  res <- read.profiling(data.dir, method = method, verbose = FALSE)

  # Probe-level taxonomies
  # (not available in the final pseq object)
  taxonomy <- res$taxonomy
  taxonomy.full <- res$taxonomy.full

  # Probe-level data
  # (not available in the final pseq object)
  probedata <- res$probedata

  res <- list(pseq = pseq, 
    	        probedata = probedata,
	        taxonomy = tax_table(as.matrix(taxonomy)),
		taxonomy.full = tax_table(as.matrix(taxonomy.full)))
  
  return(res)
  
} 

