#' @title calculate.rpa
#' @description RPA for HITChip
#' 
#' @param level level
#' @param phylo phylo
#' @param oligo.data oligo.data
#'
#' @return RPA preprocessed data
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
calculate.rpa <- function (level, phylo, oligo.data) {

  # List entities (e.g. species)
  phylo.list <- split(phylo, phylo[[level]])
  entities <- names(phylo.list)

  # initialize
  summarized.matrix <- array(NA, dim = c(length(entities), ncol(oligo.data)), dimnames = list(entities, colnames(oligo.data)))
  noise.list <- list() 

  for (entity in names(phylo.list)) {
    message(entity)

    # Pick expression for particular probes
    probes <- unique(phylo.list[[entity]][, "oligoID"])

    # oligo.data is already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 

    # dat: probes x samples
    if (nrow(dat) < 2) {
      vec <- as.vector(dat) # NOTE: circumvent RPA if there are no replicates 
      noise <- NA
    } else {
      res <- rpa.fit(dat)
      vec <- res$mu # probeset summary
      variances <- res$tau2
      names(variances) <- probes
    }

    noise.list[[entity]] <- variances
    summarized.matrix[entity, ] <- vec 
    # epsilon, alpha, beta, tau2.method, d.method

  }

  list(emat = summarized.matrix, tau2 = noise.list)
  
}