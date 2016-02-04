#' @title Maturity index
#' @description Calculate maturity index for a phyloseq object
#' @param x \code{\link{phyloseq-class}} object 
#' @return Maturity index
#' @examples
#'  #This is an artificial example based on readily available
#'  #adult data set whereas
#'  #the maturity index is typically calculated from and
#'  #used for babies/children:
#'   data("atlas1006")
#'   pseq <- atlas1006
#'   pseq <- subset_samples(pseq, DNA_extraction_method == "r" & time == 0)
#'   pseq <- transform_phyloseq(pseq, "relative.abundance")
#'   maturity(pseq) 
#' @export
#' @references To cite this R package, see citation('microbiome').
#' The microbiota maturity index has been adapted from the following papers:
#' Subramanian S et al. Nature 510:417-421, 2014.
#' Dogra S et al. mBio 6:e02419-14, 2015.
#' Korpela K et al. Nat. Comm. 7:10410, 2016.
#' @details Microbiota maturity index has been shown to differentiate healthy children (see the references). In Korpela et al. (2016) this was calculated as the first principal coordinate from a PCoA (MDS) using only significantly age-associated genus-level taxa (all groups included). In this function NMDS is used instead. The maturity index is also adjusted for age.
#' @importFrom phyloseq get_variable
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
maturity <- function (x) {

  # Pick significant age-associated taxa
  age.taxa <- rownames(lm_phyloseq(x, "age"))

  # Taxa x samples matrix for age-associated taxa
  #otu <- get_sample(x)[age.taxa,]
  pf <- prune_taxa(age.taxa, x)

  # NMDS ordination
  ord <- ordinate(pf, "NMDS", "bray")

  # Maturity index
  index <- scores(ord)[,1]

  # Adjust for age
  ages <- get_variable(x, "age")
  res <- lm(index ~ age)

  # Residuals provide the adjusted index
  index.adjusted <- residuals(res)

  index.adjusted

}

