



#' read.profiling.old
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
#' @examples 
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   dat <- read.profiling('frpa', data.dir = data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling.old <- function(method = "frpa", data.dir) {

  results <- list()

  message(paste("Reading Chip data from", data.dir))

  for (level in c("L1", "L2", "species")) {

    if (method == "frpa" && length(grep(method, dir(data.dir))) == 0) {
      warning("frpa method not available; using rpa instead")
      method <- "rpa"
    }

    f <- paste(data.dir, "/", level, "-", method, ".tab", sep = "")

    if (file.exists(f)) {
      tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
      colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]

      results[[level]] <- tab
    }
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
  tab <- polish.tax.table(tab)
  results[["taxonomy"]] <- tab  

  # Metadata      
  # Read simulated example metadata
  #library(gdata)
  #metadata.file <- paste(data.dir, "/metadata.xls", sep = "")
  #metadata <- read.xls(metadata.file, as.is = TRUE)
  #rownames(metadata) <- metadata$sampleID    
  f <- paste(data.dir, "/meta.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  rownames(tab) <- tab$sample
  results[["meta"]] <- tab      

  # Convert to numeric
  for (level in c("oligo", "species", "L1", "L2")) {

    tab <- results[[level]]        
    rnams <- rownames(tab)
    cnams <- colnames(tab)

    tab <- apply(tab, 2, as.numeric)
    rownames(tab) <- rnams
    colnames(tab) <- cnams

    if (any(is.na(tab))) {
        warning(paste("The", level, " matrix has ", sum(is.na(tab)), 
                      " missing values \n
                      - imputing.."))

        for (i in 1:ncol(tab)) {
	  inds <- which(is.na(tab[,i]))
	  if (length(inds) > 0) {
            tab[inds, i] <- sample(tab[-inds, i], length(inds))
  	  }
	}

    }

    results[[level]] <- tab        

  }
    
  results
    
} 

