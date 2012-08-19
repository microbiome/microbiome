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


#' Description: Computes the given lme model for each variable (row) in the given data frame
#'      Designed for HITchip matrices, but is applicable to any other matrix also.
#'      NOTE: does not take log automatically!
#'
#' Arguments:
#'   @param fixed TBA
#'   @param vars TBA
#'   @param d.matrix TBA
#'   @param d.info TBA
#'   @param random TBA
#'   @param correlation TBA
#'   @param weights TBA
#'   @param subset TBA
#'   @param method REML or ML
#'   @param na.action TBA
#'   @param control TBA
#'   @param contrasts TBA
#'   @param keep.data TBA
#'
#' Returns:
#'   @return lmeMatrixRes
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

lmeMatrix <- function(fixed, vars, d.matrix, d.info, random,
                        correlation = NULL, weights = NULL, subset, 
			method = c("REML", "ML"),                                                                      
                        na.action = na.fail, control = list(), contrasts = NULL,
                        keep.data = TRUE) {

  ## initialize
  vars <- rownames(d.matrix)
  lmeMatrixRes <- list()

  message("Computing lm for the following variables and data:\n")
  print(vars)
  print(str(d.info))
  print(str(t(d.matrix)))
  
  for(v in vars){
    print(v)
    data <- d.info
    data$y <- t(d.matrix[v,])

    lmeMatrixRes$lme[[v]] <- try(lme(fixed, data, random, correlation, weights, subset,
                         method, na.action, control, contrasts, keep.data))
    lmeMatrixRes$anova[[v]] <- try(anova(lmeMatrixRes$lme[[v]]))
    lmeMatrixRes$summary[[v]] <- try(summary(lmeMatrixRes$lme[[v]]))
    tmp <- try(intervals(lmeMatrixRes$lme[[v]]))

    lmeMatrixRes$intervals[[v]] <- as.list(tmp)
    
  }
  lmeMatrixRes
}




#' Description: Get Pvals From Lm Matrix
#'
#' Arguments:
#'   @param lmMatrix TBA
#'
#' Returns:
#'   @return anova.pvals
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

getPvalsFromLmMatrix <- function(lmMatrix){
  getP <- function(x){
    if(!is.element("p-value",names(x))) {
      tmp <- t(x['Pr(>F)'])
      tmp <- tmp[1:(length(tmp)-1)]
      names(tmp) <- rownames(x)[1:length(tmp)]
      return(tmp)
    } else {
      if(is.element("p-value",names(x))) {
        tmp <- t(x['p-value'])
        tmp <- tmp[2:length(tmp)]
        names(tmp) <- rownames(x)[2:(length(tmp)+1)]
        return(tmp)
      } else {
        print("Class of supposed lm/lme-object not known!")
      }
    }
  }
  anova.pvals <- sapply(lmMatrix[["anova"]], getP)
  ##rownames(anova.pvals) <- rownames(lmMatrix[[2]][[1]]['Pr(>F)'])
  return(anova.pvals)
}


