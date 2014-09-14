# Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' diversity.table
#'
#' Description: Estimate diversity within each taxonomic group across the samples
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale; this should present the highest level phylogeny (typically probe-level) used in the analysis
#'   @param phylogeny.info mapping table between taxonomic levels
#'   @param level.from higher-level taxonomic groups for which the diversity will be calculated
#'   @param level.to lower taxonomic level used for the diversity calculations. Must correspond to the input data matrix argument 'dat'
#'   @param diversity.index diversity index (shannon or invsimpson) 
#'   @param det.th Optional detection threshold. 
#'   @param min.probes Minimum number of probes for a phylotype
#' 
#' Returns:
#'   @return Table with various richness, evenness, and diversity indicators
#'
#' @examples \dontrun{divtab <- diversity.table(dat, 
#'   	     	       		level.from = "L2", 
#'				level.to = "oligo")}
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

diversity.table <- function (dat, level.from, level.to, phylogeny.info = NULL, diversity.index = "shannon", det.th = 0, min.probes = 0) {

  #phylogeny.info = phylogeny.info; level.from = "L1"; level.to = "oligo"; diversity.index = "shannon"; det.th = 0; min.probes = 0

  if (is.null(phylogeny.info)) {
    phylogeny.info <- GetPhylogeny("HITChip", "filtered")
  }

  # Ensure that the same phylogeny version used in the data and phylogeny
  # by taking the common part only
  if (level.to == "oligo") {
    coms <- intersect(phylogeny.info$oligoID, rownames(dat))
    phylogeny.info <- phylogeny.info[phylogeny.info$oligoID %in% coms,]
    dat <- dat[coms,]
  }

  if (level.to == "oligo") {level.to <- "oligoID"}		
  if (level.to == "level 1") {level.to <- "L1"}		
  if (level.to == "level 2") {level.to <- "L2"}		

  if (level.from == "level 1") {level.from <- "L1"}		
  if (level.from == "level 2") {level.from <- "L2"}		

  if (!any(rownames(dat) %in% phylogeny.info[[level.to]])) {
    stop("Provide input data matrix that corresponds to the target level ie. level.to argument!")
  }

  level.data <- levelmap(level.from = level.from, level.to = level.to, phylogeny.info = phylogeny.info)

  # Include only phylotypes with at leas min.probes
  pts <- names(which(sapply(level.data, length) >= min.probes))
  level.data <- level.data[pts]
		
  tab <- array(NA, dim = c(length(level.data), ncol = ncol(dat)))		
  rownames(tab) <- names(level.data)
  colnames(tab) <- colnames(dat)

  for (nam in names(level.data)) {

    o <- level.data[[nam]]

    if (length(o) > 1) {
      divs <- estimate.diversity(dat[o, ], diversity.index = diversity.index, det.th = det.th)
    } else {
      warning(paste("Not enough oligos for ", nam, ": diversity calculations skipped!"))
      divs <- NULL
    }    


    if (diversity.index %in% c("shannon", "invsimpson")) {
       divs <- divs$diversity
    } else if (diversity.index %in% c("richness")) {
       divs <- divs$richness
    } else if (diversity.index %in% c("evenness")) {
       divs <- divs$evenness
    }

    if (!is.null(divs)) {
      tab[nam, ] <- divs
    } else {
      tab[nam, ] <- rep(NA, ncol(tab))    
    }

  }

  tab

}


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
#'
#' @examples data(peerj32); 
#' 	     div <- estimate.diversity(10^t(peerj32$microbes), 
#'	     	    			det.th = 0)
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate.diversity <- function (dat, diversity.index = "shannon", det.th = NULL) {

  # Species diversity
  # Always use the complete data for diversity calculations
  # If you wish calculate diversity for thresholded data
  # this has to be done manually
  H <- microbiome::diversity(dat, diversity.index = diversity.index, det.th = 0)

  # richness - count phylotypes that exceed detection threshold
  # Use automated threshold determination if det.th is NULL
  S <- microbiome::richness(dat, det.th = det.th)

  # evenness - use phylotypes that exceed detection threshold
  # Use automated threshold determination if det.th is NULL
  J <- microbiome::evenness(dat, det.th = det.th)
  
  names(J) <- names(S) <- names(H) <- colnames(dat)

  data.frame(list(evenness = J, richness = S, diversity = H))
	
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
#' @examples data(peerj32); 
#' 	     div <- diversity(10^t(peerj32$microbes), 
#'	     	              det.th = 0)
#'
#' @export
#' @import vegan
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

diversity <- function (dat, diversity.index = "shannon", det.th = 0) {

  # Impute missing values
  dat <- 10^t(impute(t(log10(dat))))
  
  # Apply detection threshold
  dat.th <- dat - det.th
  dat.th[dat.th < 0] <- 0

  # Relative abundancies
  x <- relative.abundance(dat.th)

  # Calculate diversity
  if (diversity.index == "shannon") 
    x <- -x * log(x, base = exp(1))
  else x <- x * x
  if (length(dim(x)) > 1) 
    H <- apply(x, MARGIN = 2, sum, na.rm = TRUE)
  else H <- sum(x, na.rm = TRUE)
  if (diversity.index == "invsimpson") 
    H <- 1/H

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
#'
#' @examples data(peerj32); 
#' 	     rich <- richness(10^t(peerj32$microbes), 
#'	     det.th = 100)
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
#' @examples data(peerj32); 
#' 	     eve <- evenness(10^t(peerj32$microbes), 
#' 	     	    	     det.th = 100)
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


#' relative.abundance
#'
#' Description: Estimate relative abundance for each phylotype in each sample with a given threshold
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param det.th detection threshold. 
#'
#' Returns:
#'   @return Vector containing relative proportions for each phylotype in each sample 
#'
#' @examples data(peerj32); 
#' 	     relab <- relative.abundance(10^t(peerj32$microbes), 
#'	     	      		det.th = 0)
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

relative.abundance <- function (dat, det.th = 0) {

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
  dat <- apply(dat.th, 2, function (x) {x/sum(x)})

  dat
	
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
#' @examples data(peerj32); 
#' 	     abtab <- make.abundancy.table(10^t(peerj32$microbes), 
#'	     	      			det.th = 0)
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
#'   @param ylim y-axis limits
#'
#' Returns:
#'   @return Sample group list corresponding to the boxplot groups.
#' @examples data(peerj32); 
#'           div <- diversity.boxplot(peerj32$microbes, 
#'	     	                      sample.groups = list(1:22, 23:44), 
#'				      det.th = 0)
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

diversity.boxplot <- function (dat, sample.groups, diversity.index = "shannon", det.th = NULL, title = NULL, col.list = NULL, ylim = NULL) {

  div.table <- estimate.diversity(dat, diversity.index = diversity.index, det.th = det.th)

  sample.names <- colnames(dat)

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
  names(div) <- sample.names

  if (is.null(col.list)) {
    # Define colors for each sample group
    col.list <- gray(seq(0, 1, length = length(sample.groups)))
    names(col.list) <- names(sample.groups)
  }

  # Boxplot
  if (is.null(ylim)) { ylim <- c(min(div), max(div)) }
  boxplot(div[sample.groups[[1]]], col=col.list[[1]],xlim=c(0,length(sample.groups))+0.5, main = title, ylab = ylab, las = 2, ylim = ylim)
  axis(1,labels=names(sample.groups), at=c(1:length(sample.groups)))
  for (i in 2:length(sample.groups)) {
    boxplot(div[sample.groups[[i]]], col=col.list[[i]], names=names(sample.groups)[i], add = T, at = i, ylab = NULL, yaxt = "n")
  }

  list(sample.groups = sample.groups, diversity.table = div.table)

}




