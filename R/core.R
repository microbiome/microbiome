#' @title Core Microbiota
#' @description Filter the phyloseq object to include only prevalent taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold for absence/presence (greater or equal).
#' @param prevalence Prevalence threshold (in [0, 1]; greater or equal)
#' @param method Either "standard" or "bootstrap". The standard methods selects
#' the taxa that exceed the given detection and prevalence threshold.
#' The bootstrap method is more robust an described in Salonen et al. (2012).
#' Note that the results may depend on the random seed unless a sufficiently
#' large bootstrap sample size is used.
#' @param Nsample Only needed for method "bootstrap". Bootstrap sample size, default is the same size as data.
#' @param bs.iter Only needed for method "bootstrap". Bootstrap iterations.
#' @param I.max Only needed for method "bootstrap". Upper limit for intensity threshold. Later addition. Set to NULL (default) to replicate Salonen et al.
#' @param include.lowest Include the lower boundary of the detection and prevalence cutoffs. FALSE by default.
#' @param ... Arguments to pass.
#' @return Filtered phyloseq object including only prevalent taxa
#' @references
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @aliases filter_prevalent
#' @examples
#'   data(dietswap)
#'   # In practice, use more bootstrap iterations
#'   pseq <- core(dietswap, 200, .2, bs.iter = 20)
#'
core <- function (x, detection, prevalence, method = "standard", Nsample = NULL, bs.iter = 1000, I.max = NULL, include.lowest = FALSE) {

  # TODO: add optional renormalization such that the core member
  # abundances would sum up to 1 ?
  taxa <- core_members(x, detection, prevalence, method, Nsample = Nsample, bs.iter = bs.iter, I.max = I.max, include.lowest = include.lowest)
  prune_taxa(taxa, x)

}






