#' read.profiling
#' 
#' Read run.profiling.script output into R
#'
#' @param method ('frpa' / 'rpa' / 'sum' / 'ave')
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
#' @param impute impute missing oligo signals
#' 
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples 
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   dat <- read.profiling('frpa', data.dir = data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling <- function(method = "frpa", data.dir, impute = TRUE) {
    
  results <- list()

  message(paste("Reading Chip data from", data.dir))

  for (level in c("L0", "L1", "L2", "species")) {

    if (method == "frpa" && length(grep(method, dir(data.dir))) == 0) {
      warning("frpa method not available; using rpa instead")
      method <- "rpa"
    }

    f <- paste(data.dir, "/", level, "-", method, ".tab", sep = "")

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
    colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]

    results[[level]] <- tab

  }

  # oligo data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  results[["oligo"]] <- tab

  # Phylogeny
  f <- paste(data.dir, "/phylogeny.full.tab", sep = "")
  #f <- paste(data.dir, "/phylogeny.filtered.tab", sep = "")        
  tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  tab <- polish.phylogeny.info(tab)
  results[["taxonomy"]] <- tab  

  # Metadata      
  # Read simulated example metadata
  #library(gdata)
  #metadata.file <- paste(data.directory, "/metadata.xls", sep = "")
  #metadata <- read.xls(metadata.file, as.is = TRUE)
  #rownames(metadata) <- metadata$sampleID    
  f <- paste(data.directory, "/metadata.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  rownames(tab) <- tab$sample
  return(tab)
  results[["meta"]] <- tab      

  # Convert to numeric
  for (level %in% c("oligo", "species", "L0", "L1", "L2")) {
    tab <- results[[level]]        
    rnams <- rownames(tab)
    cnams <- colnames(tab)
    tab <- apply(tab, 2, as.numeric)
    rownames(tab) <- rnams
    colnames(tab) <- cnams
    
    if (impute && any(is.na(tab))) {
        warning(paste("The", level, " matrix has ", sum(is.na(tab)), 
                      " missing values \n
                      - imputing.."))
        tab <- 10^t(impute(t(log10(tab))))
    }

    results[[level]] <- tab        

  }
    
  results
    
} 

