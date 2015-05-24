#' Diversity estimation. Augments the estimate_richness function of the
#' phyloseq package
#'
#' @param x \code{\link{phyloseq-class}} object 
#' @param split See help(phyloseq::estimate_richness)
#' @param measures See help(phyloseq::estimate_richness)
#' @param det.th detection threshold for observing taxa (absence / presence)
#'
#' @return Vector containing relative proportions for each phylotype in 
#'         each sample 
#'
#' @examples 
#'   #pseq <- download_microbiome("peerj32")$physeq
#'   #d <- estimate_diversity(pseq, det.th = 0)
#'
#' @importFrom phyloseq estimate_richness
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
estimate_diversity <- function(x, split = TRUE, measures = NULL, det.th = 0) {

  res <- NULL
  measures1 <- setdiff(measures, "Richness")
  if (is.null(measures1) || length(measures1) > 0) {   
    res <- estimate_richness(x, split = split, measures = measures1)
  } 

  if (("Observed" %in% measures) || is.null(measures)) {

    # Pick the OTU data
    otu <- otu_table(x)@.Data

    # Calculate richness.
    # This simply indicates how many taxa are present in each sample
    # (exceed the detection threshold). This measure is sometimes used with
    # phylogenetic microarrays.
    ri <- colSums(otu > det.th)

    # Add to result data.frame
    if (is.null(res)) {
      res <- data.frame(Observed = ri)
    } else {
      res$Observed <- ri
    }
  }

  res

}


#' relative.abundance
#'
#' Estimate relative abundance for each phylotype in each sample 
#' with a given threshold '
#'
#' @param dat data matrix (phylotypes x samples) in original (non-log) scale
#' @param det.th detection threshold. 
#'
#' @return Vector containing relative proportions for each phylotype in 
#'         each sample 
#'
#' @examples 
#'   data(peerj32)
#'   relab <- relative.abundance(10^t(peerj32$microbes), det.th = 0)
#'
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

relative.abundance <- function(dat, det.th = 0) {
    
    # Specify detection threshold if not provided Use the 80% quantile
    # as this has proven robust across methodologies

    if (is.null(det.th)) {
        det.th <- quantile(na.omit(dat), 0.8)
        warning(paste("Applying detection threshold at 0.8 quantile: ", det.th))
    }
    
    # Apply detection threshold
    dat.th <- dat - det.th
    dat.th[dat.th < 0] <- 0
    
    # Species richness - count phylotypes that exceed detection threshold
    dat <- apply(dat.th, 2, function(x) {x <- na.omit(x); x/sum(x)})
    
    dat
    
}

#' make.abundancy.table
#'
#' Calculate abundancies
#' Discretize Hitchip matrix to form abundancy table
#' of form j, nj where j is number of counts and nj is number
#' of phylotypes with the corresponding counts
#' this format is often required by richness estimation
#'
#' @param dat data matrix
#' @param det.th detection threshold
#' @param discretization.resolution discretization resolution
#'
#' @return abundancy table
#'
#' @examples 
#'   data(peerj32)
#'   abtab <- make.abundancy.table(10^t(peerj32$microbes), det.th = 0)
#' @export
#' @references 
#'    To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

make.abundancy.table <- function(dat, det.th, discretization.resolution = 1) {
    
    di <- 10^dat - 10^det.th
    di[di < 0] <- 0
    di <- discretization.resolution * (di)
    di[di > 0 & di < 1] <- 1
    di <- round(di)
    
    ab <- t(di)
    
    ab
}

