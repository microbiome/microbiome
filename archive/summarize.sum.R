#' Probeset summarization with SUM
#' 
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param level taxonomic level for the summarization. 
#'   @param probedata preprocessed probes x samples data matrix in absolute domain
#'   @param verbose print intermediate messages
#'   @param downweight.ambiguous.probes Downweight probes with multiple targets
#'
#' Returns:
#'   @return List with two elements: abundance.table (summarized data matrix in absolute scale) and probe.parameters used in the calculations
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
summarize.sum <- function (taxonomy, level, probedata, verbose = TRUE, downweight.ambiguous.probes = TRUE) {

  # Convert to log10 domain	      
  oligo.data <- probedata
  probe.parameters <- list()
 
  probesets <- retrieve.probesets(taxonomy, level = level)

  if (downweight.ambiguous.probes) {
    nPhylotypesPerOligo <- n.phylotypes.per.oligo(taxonomy, level) 
    probe.weights <- 1/nPhylotypesPerOligo
  } else {
    probe.weights <- rep(1, nrow(taxonomy))
    names(probe.weights) <- rownames(taxonomy)
  }

  # initialize
  summarized.matrix <- matrix(NA, nrow = length(probesets), 
  		       		  ncol = ncol(oligo.data))
  rownames(summarized.matrix) <- names(probesets)
  colnames(summarized.matrix) <- colnames(oligo.data)

  for (set in names(probesets)) {

    # print(set)

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- as.matrix(oligo.data[probes,])
    if (length(probes) == 1)  {
      dat <- as.matrix(oligo.data[probes,], nrow = length(probes))
    }
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    # Weight each probe by the inverse of the number of matching phylotypes
    # Then calculate sum -> less specific probes are downweighted
    # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
    if (nrow(dat) > 1) {
      dat <- dat * probe.weights[rownames(dat)]
      vec <- colSums(dat, na.rm = T)               
    } else {
      vec <- as.vector(unlist(dat))
    }

    summarized.matrix[set, ] <- vec

  }

  list(abundance.table = summarized.matrix, probe.parameters = probe.weights)
  
}

