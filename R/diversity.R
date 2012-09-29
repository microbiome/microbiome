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
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param diversity.index diversity index (shannon or invsimpson) 
#'   @param det.th detection threshold. Used for richness and evenness estimation. Not used in diversity estimation. 
#'
#' Returns:
#'   @return Table with various richness, evenness, and diversity indicators
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate.diversity <- function (dat, diversity.index = "shannon", det.th = NULL) {

  dat.orig <- dat

  veganT <- require(vegan)
  if(!veganT) { install.packages("vegan") }

  # Species diversity
  # Always use the complete data for diversity calculations
  # If you wish calculate diversity for thresholded data
  # this has to be done manually
  H <- microbiome::diversity(dat.orig, diversity.index = diversity.index, det.th = 0)

  # richness - count phylotypes that exceed detection threshold
  # Use automated threshold determination if det.th is NULL
  S <- microbiome::richness(dat.orig, det.th = det.th)

  # evenness - use phylotypes that exceed detection threshold
  # Use automated threshold determination if det.th is NULL
  J <- microbiome::evenness(dat.orig, det.th = det.th)
  
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


#' diversity
#'
#' Description: Estimate diversity for each sample with a given threshold
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param diversity.index diversity index (shannon or invsimpson) 
#'   @param det.th detection threshold. Used for richness and evenness estimation. Not used in diversity estimation. 
#'
#' Returns:
#'   @return Vector containing diversity estimate for each sample 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

diversity <- function (dat, diversity.index = "shannon", det.th = 0) {

  veganT <- require(vegan)
  if(!veganT) { install.packages("vegan") }

  # Specify detection threshold if not provided
  # Use the 80% quantile as this has proven robust across methodologies
  if (is.null(det.th)) {
    dat <- 10^t(impute(t(log10(dat)))) # impute missing values
    det.th <- quantile(dat, 0.8)
    warning(paste("Applying detection threshold at 0.8 quantile: ", det.th))
  }

  # Apply detection threshold
  dat.th <- dat - det.th
  dat.th[dat.th < 0] <- 0

  # Species diversity
  H <- vegan::diversity(dat.th, index = diversity.index, MARGIN = 2)
  names(H) <- colnames(dat)

  H
	
}


#' richness
#'
#' Description: Estimate richness for each sample with a given threshold
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param det.th detection threshold. Used for richness and evenness estimation. Not used in diversity estimation. 
#'
#' Returns:
#'   @return Vector containing richness estimate for each sample 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

richness <- function (dat, det.th = NULL) {

  # Specify detection threshold if not provided
  # Use the 80% quantile as this has proven robust across methodologies
  if (is.null(det.th)) {
    dat <- 10^t(impute(t(log10(dat)))) # impute missing values
    det.th <- quantile(dat, 0.8)
    warning(paste("Applying detection threshold at 0.8 quantile: ", det.th))
  }

  # Apply detection threshold
  dat.th <- dat - det.th
  dat.th[dat.th < 0] <- 0

  # Species richness - count phylotypes that exceed detection threshold
  S <- colSums(dat.th > 0)
  names(S) <- colnames(dat)

  S
	
}



#' evenness
#'
#' Description: Estimate evenness for each sample with a given threshold
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param det.th detection threshold. Used for richness and evenness estimation. Not used in diversity estimation. 
#'
#' Returns:
#'   @return Vector containing evenness estimate for each sample 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

evenness <- function (dat, det.th = NULL) {

  # Specify detection threshold if not provided
  # Use the 80% quantile as this has proven robust across methodologies
  if (is.null(det.th)) {
    dat <- 10^t(impute(t(log10(dat)))) # impute missing values
    det.th <- quantile(dat, 0.8)
    warning(paste("Applying detection threshold at 0.8 quantile: ", det.th))
  }

  # Apply detection threshold
  dat.th <- dat - det.th
  dat.th[dat.th < 0] <- 0

  # Species richness - count phylotypes that exceed detection threshold
  S <- colSums(dat.th > 0)
  
  # Pielou's evenness (J); always calculated from H.shannon which has to use same detection threshold
  # as the richness measure:
  H.shannon.th <- microbiome::diversity(dat.th, diversity.index = "shannon", det.th = 0)

  J <- H.shannon.th/log(S)
  names(J) <- colnames(dat)

  J
	
}



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
#'   @param dat data matrix in the original (non-log) scale: features x samples.
#'   @param sample.groups specify the distinct sample groups     
#'   @param diversity.index Options: "shannon" / "invsimpson" / "richness" / "evenness" 
#'   @param det.th Detection threshold for richness and evenness calculations
#'   @param title figure title
#'   @param col.list Optional: colors for the sample groups
#'
#' Returns:
#'   @return Sample group list corresponding to the boxplot groups.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

diversity.boxplot <- function (dat, sample.groups, diversity.index = "shannon", det.th = NULL, title = NULL, col.list = NULL) {

  div.table <- estimate.diversity(dat, diversity.index = diversity.index, det.th = det.th)

  if (is.null(title)) {title <- diversity.index}

  if (diversity.index %in% c("shannon", "invsimpson", "diversity")) {
    div <- div.table[["diversity"]]
    ylab = "Diversity"
  } else if (diversity.index %in% c("richness")) {
    div <- div.table[["richness"]]
    ylab <- "Richness"
  } else if (diversity.index %in% c("evenness")) {
    div <- div.table[["evenness"]]
    ylab <- "Evenness"
  }
  names(div) <- rownames(div.table)

  if (is.null(col.list)) {
    # Define colors for each sample group
    col.list <- gray(seq(0, 1, length = length(sample.groups)))
    names(col.list) <- names(sample.groups)
  }

  # Boxplot
  boxplot(div[sample.groups[[1]]], col=col.list[[1]], ylim=c(min(div),max(div)),xlim=c(0,length(sample.groups))+0.5, main = title, ylab = ylab, las = 2)
  axis(1,labels=names(sample.groups), at=c(1:length(sample.groups)))
  for (i in 2:length(sample.groups)) {
    boxplot(div[sample.groups[[i]]], col=col.list[[i]], names=names(sample.groups)[i], add = T, at = i, ylab = NULL, yaxt = "n")
  }

  list(sample.groups = sample.groups, diversity.table = div.table)

}




