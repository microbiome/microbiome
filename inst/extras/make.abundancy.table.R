#' @title Abundance table
#' @description Calculate abundancies.
#' Discretize Hitchip matrix to form abundancy table
#' of form j, nj where j is number of counts and nj is number
#' of phylotypes with the corresponding counts
#' this format is often required by richness estimation
#' @param dat data matrix
#' @param det.th detection threshold
#' @param discretization.resolution discretization resolution
#' @return abundancy table
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

