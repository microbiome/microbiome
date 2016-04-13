#' @title Relative abundance
#' @description Estimate relative abundance for each phylotype in each sample 
#' with a given threshold.
#' @param dat data matrix (phylotypes x samples) in original (non-log) scale
#' @param det.th detection threshold. 
#' @return Vector containing relative proportions for each phylotype in 
#'         each sample 
#' @examples 
#'   data(peerj32)
#'   relab <- relative.abundance(10^t(peerj32$microbes), det.th = 0)
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
