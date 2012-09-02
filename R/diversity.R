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
  # Use the 80% quantile as this has proven robust across methodologies
  if (is.null(det.th)) {
    dat <- 10^t(impute(t(log10(dat)))) # impute missing values
    #det.th <- 10^estimate.min.threshold(log10(dat))
    det.th <- quantile(dat, 0.8)
    warning(paste("Applying detection threshold at 0.8 quantile: ", det.th))
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



#' diversity.boxplot
#'
#' Description: diversity.boxplot
#'
#' Arguments:
#'   @param div.table diversity table: output from the estimate.diversity function
#'   @param diversity.index specify the diversity index to plot (diversity / richness / evenness)
#'   @param title figure title
#'   @param n.groups max number of groups for boxplot
#'
#' Returns:
#'   @return Sample group list corresponding to the boxplot groups.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

diversity.boxplot <- function (div.table, diversity.index = "invsimpson", title = NULL, col.list, sample.groups) { 

  if (diversity.index %in% c("shannon", "invsimpson", "diversity")) {
    div <- div.table[["diversity"]]
  } else if (diversity.index %in% c("richness")) {
    div <- div.table[["richness"]]
  } else if (diversity.index %in% c("evenness")) {
    div <- div.table[["evenness"]]
  }
  names(div) <- rownames(div.table)

  ## "Box-plotting":
  # Diversity index
  boxplot(div[sample.groups[[1]]], ylab = diversity.index, col=col.list[[1]], ylim=c(min(div),max(div)),xlim=c(0,length(sample.groups))+0.5, main = title)
  axis(1,labels=names(sample.groups),at=c(1:length(sample.groups)))
  for (i in 2:length(sample.groups)) {
	boxplot(div[sample.groups[[i]]],col=col.list[[i]],names=names(sample.groups)[i],add=T,at=i)
  }

  sample.groups

}




