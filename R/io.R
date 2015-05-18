#' read.profiling
#' 
#' Read run.profiling.script output into R
#'
#' @param method ('frpa' / 'rpa' / 'sum' / 'ave')
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
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

read.profiling <- function(data.dir) {

  results <- list()

  message(paste("Reading Chip data from", data.dir))

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  otu <- t(tab)
  #results[["oligo"]] <- tab

  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  # Convert into phyloseq taxonomyTable format
  taxonomy <- tax_table(as.matrix(tab))     
  #results[["taxonomy"]] <- tab

  # Read sample metadata      
  f <- paste(data.dir, "/meta.tab", sep = "")
  if (file.exists(f)) {
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    rownames(tab) <- tab$sample
    meta <- tab
    #results[["meta"]] <- tab      
  }

  # Convert the object into phyloseq format
  res <- hitchip2physeq(otu, meta, taxonomy, detection.limit = 0)

  res
    
} 


