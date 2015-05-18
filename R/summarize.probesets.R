#' Description: summarize.probesets
#'
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param level summarization level
#'   @param verbose print intermediate messages
#'   @param species.matrix Optional. Provide pre-calculated species-level summaries to speed up computation.
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
summarize.probesets <- function (taxonomy, oligo.data, method, level, verbose = TRUE, species.matrix = NULL) {

  if (level == "species") {
    method <- gsub(".through.species", "", method)
  } 

  if (method %in% c("rpa", "frpa", "rpa.with.affinities", "sum.through.species", "ave.through.species")) {

    # Summarize probes through species level (default with RPA)
    res <- summarize.probesets.through.species(level = level, taxonomy = taxonomy, oligo.data = oligo.data, method = gsub(".through.species", "", method), verbose = verbose)

    } else if (method %in% c("rpa.direct", "rpa.with.affinities.direct", "sum", "ave")) {

    # Option 2: Summarize from oligos to all levels directly (default with SUM and AVE)
    res <- summarize.probesets.directly(level, taxonomy, oligo.data, gsub(".direct", "", method))

  } else {
    stop(paste("method", method, "not implemented in HITChipDB"))
  }
  
  res

}

