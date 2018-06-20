#' @title Read profiling
#' @description Read run.profiling.script output into R
#' @param level phylogenetic level ('oligo' / 'species' / 'L1' / 'L2' / 'L0') 
#' 	  	or 'phylogeny.full', 'phylogeny.filtered'
#' @param method ('frpa' / 'rpa' / 'sum' / 'ave')
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
#' @param log10 Logical. Logarithmize the data TRUE/FALSE. 
#'              By default, the data is in original non-log scale.
#' @param impute impute missing oligo signals
#' @return data matrix (phylo x samples)
#' @export
#' @examples
#'  \dontrun{
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   dat <- read.profiling('L1', 'frpa', data.dir = data.dir)
#'  }
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling <- function(level = NULL, method = "frpa", data.dir, 
                         log10 = TRUE, impute = TRUE) {

    .Deprecated("read_hitchip")
    
    # level <- 'oligo'; method = 'sum'; data.dir = 'test/'; log10 = TRUE
    if (level %in% c("L0", "L1", "L2", "species")) {
        if (method == "frpa" && length(grep(method, dir(data.dir))) == 0) {
            warning("frpa method not available; using rpa instead")
            method <- "rpa"
        }
        f <- paste(data.dir, "/", level, "-", method, ".tab", sep = "")

    } else if (level == "oligo") {
        f <- paste(data.dir, "/oligoprofile.tab", sep = "")
    } else if (level == "phylogeny.full") {
        f <- paste(data.dir, "/phylogeny.full.tab", sep = "")
    } else if (level %in% c("phylogeny.filtered")) {
        f <- paste(data.dir, "/phylogeny.filtered.tab", sep = "")
    } else if (level %in% c("phylogeny.info", "phylogeny.full", "phylogeny")) {
        f <- paste(data.dir, "/phylogeny.full.tab", sep = "")
    }
    
    
    message(paste("Reading", f))
    
    if (level %in% c("L0", "L1", "L2", "species")) {
        
        tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, 
                        as.is = TRUE)
        colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
        
    } else if (level == "oligo") {
        
        tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, 
                        as.is = TRUE)
        colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
        
    } else if (length(grep("phylogeny", level)) > 0) {
        
        tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
        
    }
    
    # Convert to numeric
    if (level %in% c("oligo", "species", "L0", "L1", "L2")) {
        
        rnams <- rownames(tab)
        cnams <- colnames(tab)
        tab <- apply(tab, 2, as.numeric)
        rownames(tab) <- rnams
        colnames(tab) <- cnams
        
    }
    
    if (impute && any(is.na(tab))) {
        warning(paste("The matrix has ", sum(is.na(tab)), 
                      " missing values \n                             
                      - imputing."))
        tab <- 10^t(impute(t(log10(tab))))
    }
    
    if (log10 && (level %in% c("oligo", "species", "L0", "L1", "L2"))) {
        message("Logarithmizing the data")
        tab <- log10(tab)
    }
    
    tab
    
} 

