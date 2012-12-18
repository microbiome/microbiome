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


#' Description: count
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param d TBA
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

count <- function(d){
   tabulate(d)
}



#' Description: summarize.oligos 
#'
#' Arguments:
#'   @param oligo.matrix oligo.matrix
#'   @param phylogeny.info phylogeny.info
#'   @param level taxonomic level
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.oligos <- function(oligo.matrix, phylogeny.info, level = "L2"){

  phylogeny.info <- polish.phylogeny.info(phylogeny.info)
		 
  oligo <- split(phylogeny.info$oligoID,phylogeny.info[,level])
  oligo <- lapply(oligo, function(x) unique(as.character(x)))

  Sim <- sapply(oligo,function(x){  
    if (length(x) > 1)
      return(colSums(oligo.matrix[x, ]))
    else
      return(oligo.matrix[x, ])
  })
  t(Sim)
}


#' Description: mixingMatrix
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param phylogeny.info phylogeny.info
#'   @param level taxonomic level
#'
#' Returns:
#'   @return oligos x phylotypes mixing matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

mixingMatrix <- function(phylogeny.info, level){

  M <- matrix(0,length(unique(phylogeny.info$oligoID)),length(unique(phylogeny.info[,level])),dimnames=list(sort(as.character(unique(phylogeny.info$oligoID))),sort(as.character(unique(phylogeny.info[,level])))))

  for (i in 1:nrow(phylogeny.info))
    M[as.character(phylogeny.info$oligoID[i]),as.character(phylogeny.info[i,level])]=1

  M <- apply(M, 2, function(x) x/sum(x))

  return(M)

}

#' Description: ngp: Non-negative probabilistic solution to probe mixing problem.
#'
#' Uses gamma prior to solve known cross-hybridisation problems. 
#' Perfect probe targets are given in phylogeny.info data frame. 
#'
#' Arguments:
#'   @param oligo.data oligo.data in absolute domain
#'   @param phylogeny.info oligo - phylotype mapping data frame
#'   @param level taxonomic level
#'   @param lambda - stregth of gamma prior. Default: 0.001.
#'   @param alpha - alpha parameter of gamma prior. Default: 1.
#'   @param beta - beta parameter of gamma prior. Default: 1.
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

ngp <- function(oligo.data, phylogeny.info, level, lambda=0.001, alpha=1,beta=1){

   # oligos x phylotypes mixing matrix for the given level
   M <- mixingMatrix(phylogeny.info, level)
   coms <- intersect(rownames(oligo.data), rownames(M))
   M <- M[coms, ]  
   oligo.data <- oligo.data[coms, ]

   # starting guess using pseudoinverse
   require(MASS)
   W=t(M)%*%M
   X=t(oligo.data) %*% M
   A=ginv(W)%*%t(X)
   A[which(A<0)]=0

   for (i in 1:ncol(oligo.data)){
      Acol=A[,i]
      for (j in 1:ncol(M)){  
         a=2*W[j,j]
         b=sum((W[j,-j]+W[-j,j])*Acol[-j])-2*X[i,j]-beta*lambda
         c=lambda*(alpha-1)
         if (a!=0){
           D=b^2-4*a*c
           if (D>0)
             A[j,i]=max((-b+sqrt(D)),(-b-sqrt(D)))/(2*a)
         }else
            A[j,i]=-c/b
     }
   }
   colnames(A) <- colnames(oligo.data)
   rownames(A) <- colnames(M)
   return(A)
}

#' Description: Deconvolution
#'
#' !OBSOLETE! Used nmf package (currently removed from CRAN) for solving cross-hyb problem.
#' Calls function ngp with standard arguments (SEE: ngp).
#'
#' Arguments:
#'   @param oligo.data oligo.data in absolute domain
#'   @param phylogeny.info oligo - phylotype mapping data frame
#'   @param level taxonomic level
#'   @param block.solution block.solution
#'   @param verbose verbose
#'   @param ... parameters to be passed to ngp function
#'
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

deconvolution.nonneg <- function(oligo.data, phylogeny.info, level, block.solution = T,verbose=F,...){
   ngp(oligo.data, phylogeny.info, level,...)
}

# --------------------------------------------------------------------



