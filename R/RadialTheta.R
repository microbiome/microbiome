#' @title Radial theta function.
#' @description Adapted from \pkg{NeatMap} and \pkg{phyloseq} packages. In phyloseq this was vectorized for speed and simplicity, also only calculates theta and not r. Copied here as it is not exported and hence not available via phyloseq
#' @keywords internal
RadialTheta <- function(pos){
    pos <- as(pos, "matrix")
    xc  <- mean(pos[, 1])
    yc  <- mean(pos[, 2])
    theta <- atan2((pos[, 2] - yc), (pos[, 1] - xc))
    names(theta) <- rownames(pos)
    return(theta)
}