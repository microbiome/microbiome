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
#'   # dat <- read_profiling(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read.profiling <- function(data.dir, method) {

  message(paste("Reading Chip data from", data.dir))
  results <- list()

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  results[["probedata"]] <- tab

  # Read abundance tables
  for (s in c("L1", "L2", "species")) {
    f <- paste(data.dir, "/", s, "-", method, ".tab", sep = "")
    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
    colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
    results[[s]] <- tab
  }
  
  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  taxonomy <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  results[["taxonomy"]] <- taxonomy
    
  # Read unfiltered taxonomy table
  f <- paste(data.dir, "/taxonomy.full.tab", sep = "")
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


#' Read run.profiling.script output into R and preprocess into phyloseq format
#'
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
#' @param output Specify the desired output.
#' @param method Probe summarization method ("rpa" or "sum")
#' 
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples 
#'   # data.dir <- system.file("extdata", package = "microbiome")
#'   # dat <- read_profiling(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read_hitchip <- function(data.dir, output = "all", method = "rpa") {

  # Read	     
  res <- read.profiling(data.dir, method = method)
  taxonomy <- res$taxonomy
  taxonomy.full <- res$taxonomy.full
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

  pseq <- hitchip2physeq(t(abu), meta, taxonomy, detection.limit = 10^1.8)
  if (output == "phyloseq") {
    return(pseq) 
  }  

  if (output == "all") {
    res <- list(pseq = pseq, abundance.table = abu, 
    	        meta = meta,
    	        probedata = probedata,
	        taxonomy = tax_table(as.matrix(taxonomy)),
		taxonomy.full = tax_table(as.matrix(taxonomy.full)))
    return(res)
  }

  NULL
    
} 

