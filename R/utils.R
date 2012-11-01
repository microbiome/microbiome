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


#' PhylotypeRatios
#'
#' Calculate phylotype ratios (eg. Bacteroides vs. Prevotella etc.) for a given
#' phylotypes vs. samples matrix
#'
#' @param dat phylotypes vs. samples data matrix in log10 scale
#'
#' @return phylotype pairs x samples matrix indicating the ratio (in log10 domain) between each unique pair
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples #
#' @keywords utilities

PhylotypeRatios <- function (dat) {

  phylogroups <- rownames(dat)
  Nratios <- (length(phylogroups)^2 - length(phylogroups))/2
  Nsamples <- ncol(dat)
  ratios <- list()
  for (i in 1:(length(phylogroups)-1)) {
    for (j in (i+1):length(phylogroups)) {
      pt1 <- phylogroups[[i]]
      pt2 <- phylogroups[[j]]
      ratios[[paste(pt1, pt2, sep = "-")]] <- dat[pt1,] - dat[pt2,]
    } 
  }
  ratios <- do.call(cbind, ratios)

  t(ratios)
}


#' Center data matrix.
#' 
#' Center data matrix to 0 for each variable by removing the means.
#' 
#' 
#' @usage centerData(X, rm.na = TRUE, meanvalue = NULL)
#' @param X The data set: samples x features. Each feature will be centered.
#' @param rm.na Ignore NAs.
#' @param meanvalue Can be used to set a desired center value. The default is
#' 0.
#' @return Centered data matrix.
#' @note Note that the model assumes samples x features matrix, and centers
#' each feature.
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @references See citation("microbiome").
#' @keywords utilities maths
#' @export
#' @examples
#' 
#' #centerData(X)
#' 
centerData <- function (X, rm.na = TRUE, meanvalue = NULL) {

  # Shift data matrix (columns) to zero, or given 'meanvalue'
  
  if (!rm.na) {
    xcenter <- colMeans(X)
    X2 <- X - rep(xcenter, rep.int(nrow(X), ncol(X)))
  } else {	
    X2 <- array(NA, dim = c(nrow(X), ncol(X)), dimnames = dimnames(X))
    for (i in 1:ncol(X)) {
      x <- X[,i]
      nainds <- is.na(x)
      xmean <- mean(x[!nainds])
      X2[!nainds,i] <- x[!nainds] - xmean 	
    }
    dimnames(X2) <- dimnames(X)
  }

  if (!is.null(meanvalue)) {
    # Shift the data so that mean gets a specified value
    X2 <- X2 + meanvalue
  }


  
  X2
}


#' matrix.qvalue
#'
#' Calculate qvalues for a matrix of pvalues which may contain missing values.
#' 
#'
#' @param pvals p-value matrix
#'
#' @return q-value matrix
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples #qvals <- matrix.qvalue(pvals)
#' @keywords utilities

matrix.qvalue <- function (pvals) {

  pvec <- as.vector(pvals)
  nai <- is.na(pvec);
  qvec <- rep(NA, length(pvec))
  qvec[!nai] <- qvalue::qvalue(pvec[!nai], pi0.method = "bootstrap")$qvalue
  qmat <- matrix(qvec, nrow = nrow(pvals))
  dimnames(qmat) <- dimnames(pvals)
  qmat

}

#' polish.phylogeny.info
#'
#' Ensure phylogeny.info is in correct format
#' 
#' @param phylogeny.info phylogeny.info data frame
#'
#' @return polished phylogeny.info
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples #polished.phylogeny.info <- impute(phylogeny.info) 
#' @keywords utilities

polish.phylogeny.info <- function (phylogeny.info) {

  colnames(phylogeny.info)[which(colnames(phylogeny.info) == "level.0")] <- "L0"
  colnames(phylogeny.info)[which(colnames(phylogeny.info) == "level.1")] <- "L1"
  colnames(phylogeny.info)[which(colnames(phylogeny.info) == "level.2")] <- "L2"

  colnames(phylogeny.info)[which(colnames(phylogeny.info) == "level 0")] <- "L0"
  colnames(phylogeny.info)[which(colnames(phylogeny.info) == "level 1")] <- "L1"
  colnames(phylogeny.info)[which(colnames(phylogeny.info) == "level 2")] <- "L2"

  phylogeny.info

}


#' Impute missing values from a Gaussian. 
#' 
#' @param X data matrix (features x samples)
#'
#' @return imputed data matrix
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples #X2 <- impute(X) 
#' @keywords utilities

impute <- function (X) {

  for (i in 1:ncol(X)) {
    x <- X[, i]
    nas <- is.na(x)
    X[nas, i] <- rnorm(sum(nas), mean(x[!is.na(x)]), sd(x[!is.na(x)]))
  }
  
  X
  
}



#' Strip string i.e. remove spaces from the beginning and end
#' @param s string or character vector
#'
#' @return Stripped string
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples 
#' #s2 <- Strip(s) 
#' @keywords utilities
Strip <- function (s) {

  ss <- c()
  
  for (i in 1:length(s)) {

    si <- s[[i]]
    if (!is.na(si)) {
      # Strip string i.e. remove spaces from the beginning and end
      while (substr(si,1,1)==" ") {
        si <- substr(si, 2, nchar(si))
      }
      while (substr(si, nchar(si), nchar(si))==" ") {
        si <- substr(si, 1, nchar(si) - 1)
      }
    }
    ss[[i]] <- si
  }

  ss
}



#' Description: Sort data frame dd by columns like: esort(dd, -z, b)
#'
#' Arguments:
#'   @param x data frame to sort
#'   @param sortvar sorted variable/s
#'   @param ... further parameters to pass
#'
#' Returns:
#'   @return sorted data frame
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

esort <- function(x, sortvar, ...) {

  attach(x)
  x <- x[with(x,order(sortvar,...)),]
  return(x)
  detach(x)
}




#' Description: 
#' Get lower triangle of a square matrix 
#' as a numeric vector such that
#' row-by-row, picking elements in the order
#' 2,1;3,1;3,2;4,1,...
#'        
#' Arguments:
#'   @param mat data matrix
#' Returns:
#'   @return lower triangle as vector 
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
lower.triangle <- function (mat) {

  elements <- c()
  nr <- dim(mat)[[1]]
  nc <- dim(mat)[[2]]

  for (i in 2:nr) {
    for (j in 1:(i-1)) {
      elements<-c(elements,mat[i,j])
    }
  }       
  elements
}


