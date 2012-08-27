# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' estimate.diversity
#'
#' Description: Estimate diversities for each sample (column)
#' Aliases: get.diversity.estimates. Also replaces the function ST: diversity.indices
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param diversity.index diversity index
#'   @param det.th detection threshold. Used for richness and evenness estimation. Not used in diversity estimation. 
#'
#' Returns:
#'   @return diversity indicators
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate.diversity <- function (dat, diversity.index = "shannon", det.th = NULL) {

  veganT <- require(vegan)
  if(!veganT) { install.packages("vegan") }

  # Specify detection threshold if not provided
  if (is.null(det.th)) {
    dat <- impute(dat) # impute the very few (3) missing values
    det.th <- 10^estimate.min.threshold(log10(dat))
    warning(paste("Applying detection threshold: ", det.th))
  }

  # Apply detection threshold
  dat <- dat - det.th
  dat[dat < 0] <- 0

  # Species diversity - only for species that exceeded the thresholding above
  H <- diversity(dat, index = diversity.index, MARGIN = 2)
  
  # Species richness - count phylotypes that exceed detection threshold
  S <- colSums(dat > 0)
  
  # Pielou's evenness (J) 
  H.shannon <- diversity(dat, index = "shannon", MARGIN = 2)
  J <- H.shannon/log(S)
  
  names(J) <- names(S) <- names(H) <- colnames(dat)

  data.frame(list(evenness = J, richness = S, diversity = H, det.th = det.th))
	
}

  # Compare to ACE estimate with different discretization bins
  # -> no clear way to discretize, will affect the results quite much
  # -> skip so far
  #par(mfrow=c(3,3))
  #for (discretization.resolution in 10^(-seq(-2,6,length=9))) {
  #  ab <- make.abundancy.table(dat, det.th, discretization.resolution)
  #  S2 <- estimateR(ab)["S.ACE", ]
  #  plot(S, S2, main = discretization.resolution)
  #}
  # Chao and ACE estimators by modifying this: requiring counts, add later?
  # Estimater(floor(10^t(dat))) 
  #S <- unlist(mclapply(1:ncol(dat), function(k) {'  
  #ab <- make.abundancy.table(dat[, k], discretization.resolution = 1e-4)
    #ChaoLee1992(ab, t=10, method="ACE")$Nhat
  #}))



#' make.abundancy.table
#'
#' Description: Calculate abundancies
#' Discretize Hitchip matrix to form abundancy table
#' of form j, nj where j is number of counts and nj is number
#' of phylotypes with the corresponding counts
#' this format is often required by richness estimation
#'
#' Arguments:
#'   @param dat data matrix
#'   @param det.th detection threshold
#'   @param discretization.resolution discretization resolution
#'
#' Returns:
#'   @return abundancy table
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

make.abundancy.table <- function (dat, det.th, discretization.resolution = 1) {

  di <- 10^dat - 10^det.th
  di[di<0] <- 0
  di <- discretization.resolution*(di)
  di[di>0 & di<1] <- 1
  di <- round(di)

  ab <- t(di)

  ab
}


