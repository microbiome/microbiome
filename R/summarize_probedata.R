#' Summarize phylogenetic microarray probe-level data from given input folder.
#' 
#' @param data.dir Data folder.
#' @param probedata probe-level data matrix
#' @param taxonomy probe taxonomy
#' @param level Summarization level
#' @param method Summarization method
#' @param probe.parameters Precalculater probe parameters. Optional.
#'
#' @return data matrix (taxa x samples)
#'
#' @export
#' @examples #
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
summarize_probedata <- function(data.dir = NULL, probedata = NULL, taxonomy = NULL, level, method, probe.parameters = NULL) {

  # If the data is not given as input, read it from the data directory		     
  # message(paste("Reading Chip data from", data.dir))

  # Read probe-level data
  if (is.null(probedata)) {
    f <- paste(data.dir, "/oligoprofile.tab", sep = "")
    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
    colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
    probedata <- tab
  }

  # Read taxonomy table
  if (is.null(taxonomy)) {
    f <- paste(data.dir, "/taxonomy.tab", sep = "")
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    # Convert into phyloseq taxonomyTable format
    taxonomy <- tax_table(as.matrix(tab))     
  }

  # Summarize probes through species level
  # probeset.summaries <- summarize.probesets.species(taxonomy, oligo.data, method)
  if (method %in% c("rpa", "frpa")) {
    otu <- summarize.rpa(taxonomy, level, probedata, verbose = TRUE, probe.parameters = probe.parameters)$abundance.table
  } else if (method == "sum") {
    otu <- summarize.sum(taxonomy, level, probedata, verbose = TRUE, downweight.ambiguous.probes = TRUE)$abundance.table
  }

  otu
    
}

