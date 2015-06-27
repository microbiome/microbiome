#' Read HITChip run.profiling.script output into R
#'
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
#' @param method Select the preprocessing method that you like to check
#' @param verbose verbose
#' 
#' @return data matrix (phylo x samples)
#'
#' @examples 
#'   # data.dir <- system.file("extdata", package = "microbiome")
#'   # dat <- read_profiling(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read_profiling <- function(data.dir, method, verbose = TRUE) {

  if (verbose) {message(paste("Reading Chip data from", data.dir))}

  results <- list()

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  results[["probedata"]] <- tab

  # Read abundance tables
  #for (s in c("L0", "L1", "L2", "species")) {
  #  f <- paste(data.dir, "/", s, "-", method, ".tab", sep = "")
  #  if (file.exists(f)) {
  #    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  #    colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  #    results[[s]] <- tab
  #  }
  #}
  
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


