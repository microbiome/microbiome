#' Description: Probeset summarization with RPA
#' 
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param level taxonomic level for the summarization. 
#'   @param probedata preprocessed probes x samples data matrix in absolute domain
#'   @param verbose print intermediate messages
#'   @param probe.parameters Optional. If probe.parameters are given,
#'          the summarization is based on these and model parameters are not
#' 	    estimated. A list. One element for each probeset with the following probe vectors: 
#'	    affinities, variances
#' Returns:
#'   @return List with two elements: abundance.table (summarized data matrix in absolute scale) and probe.parameters (RPA probe level parameter estimates)
#'
#' @export
#' @importFrom RPA d.update.fast 
#' @importFrom RPA rpa.fit
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
summarize.rpa <- function (taxonomy, level, probedata, verbose = TRUE, probe.parameters = NULL) {

  # Convert to log10 domain	      
  oligo.data <- log10(probedata) 
  probeinfo <- list()
 
  probesets <- retrieve.probesets(taxonomy, level = level)
  # probesets <- probesets[setdiff(names(probesets), rm.species)]
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(taxonomy, level) 

  # initialize
  summarized.matrix <- matrix(NA, nrow = length(probesets), 
  		       		  ncol = ncol(oligo.data))
  rownames(summarized.matrix) <- names(probesets)
  colnames(summarized.matrix) <- colnames(oligo.data)

  for (set in names(probesets)) {

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

    if (length(probe.parameters) > 0) {

      # Summarize with pre-calculated variances
      vec <- d.update.fast(dat, probe.parameters[[set]])

    } else {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      alpha <- 1 + 0.1*ncol(dat)/2
      beta <- 1 + 0.1*ncol(dat)*nPhylotypesPerOligo[probes]^2
      res <- rpa.fit(dat, alpha = alpha, beta = beta)
      vec <- res$mu
      probeinfo[[set]] <- res$tau2

    }
      
    summarized.matrix[set, ] <- vec 

  }

  if (!is.null(probe.parameters)) {
    probeinfo <- probe.parameters
  }

  # Return the data in absolute scale					
  summarized.matrix <- 10^summarized.matrix

  list(abundance.table = summarized.matrix, probeinfo = probeinfo)
  
}

