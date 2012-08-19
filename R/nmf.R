# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package

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
#' @author Contact: Jarkko Salojarvi \email{jarkko.salojarvi@@helsinki.fi}
#' @keywords utilities

count <- function(d){
  sapply(unique(d),function(x) sum(d==x))
}


#' Description: simulate.hitchip
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param PH.i oligomap
#'   @param Ns Ns
#'   @param level phylogenetic level
#'   @param N N
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{jarkko.salojarvi@@helsinki.fi}
#' @keywords utilities

simulate.hitchip <- function(PH.i, Ns, level = "level 2", N = 5000){

  oligo <- split(PH.i$oligoID,PH.i[, level])
  oligo <- lapply(oligo,function(x) unique(as.character(x)))
  oligos <- unique(as.character(PH.i$oligoID))
  M <- matrix(0,length(oligos),Ns)
  rownames(M) <- oligos
  s.vec <- t(round(replicate(length(oligo),sapply(rnorm(Ns,mean=N,sd=N/4),function(y) max(y,1)))))
  for (i in 1:Ns){
    for (j in 1:nrow(s.vec)){
      a <- count(sample(oligo[[j]],s.vec[j,i],replace=T))
      M[names(a),i] <- M[names(a),i] + a
   } 
  }
  rownames(s.vec) <- names(oligo)
  list(Simulated = M[order(rownames(M)),], True = s.vec)
}


#' Description: summarize.oligos
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param Data Data
#'   @param PH.i PH.i
#'   @param level phylogenetic level
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{jarkko.salojarvi@@helsinki.fi}
#' @keywords utilities

summarize.oligos <- function(Data, PH.i, level = "level.2"){
  oligo <- split(PH.i$oligoID,PH.i[,level])
  oligo <- lapply(oligo,function(x) unique(as.character(x)))
  Sim <- sapply(oligo,function(x){  
    if (length(x)>1)
      return(colSums(Data[x,]))
    else
      return(Data[x,])
  })
  t(Sim)
}


#' Description: mixingMatrix
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param oligo.map oligo.map
#'   @param level phylogenetic level
#'
#' Returns:
#'   @return oligos x phylotypes mixing matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{jarkko.salojarvi@@helsinki.fi}
#' @keywords utilities

mixingMatrix <- function(oligo.map,level){

  M <- matrix(0,length(unique(oligo.map$oligoID)),length(unique(oligo.map[,level])),dimnames=list(sort(as.character(unique(oligo.map$oligoID))),sort(as.character(unique(oligo.map[,level])))))

  for (i in 1:nrow(oligo.map))
    M[as.character(oligo.map$oligoID[i]),as.character(oligo.map[i,level])]=1

  M <- apply(M, 2, function(x) x/sum(x))

  return(M)

}



#' Description: mixingMatrix
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param oligo.data oligo.data
#'   @param oligo.map oligo.map
#'   @param level phylogenetic level
#'   @param block.solution block.solution
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{jarkko.salojarvi@@helsinki.fi}
#' @keywords utilities
deconvolution.nonneg <- function(oligo.data, oligo.map, level, block.solution = T){

   require(NMF)

   # oligos x phylotypes mixing matrix for the given level
   M <- mixingMatrix(oligo.map, level)
   M <- M[rownames(oligo.data), ]

   #least squares solution, may contain negative values.
   #H=solve(t(M) %*% M, t(M) %*% dd$Simulated)
   #this is a faster version. (but not much)
   #H.nonneg=NMF:::.fcnnls(t(M) %*% M,t(M) %*% oligo.data,pseudo=TRUE)

   if (block.solution){

     require(ggm)
     A <- t(M) %*% M
     B <- t(M) %*% oligo.data
     H.nonneg <- matrix(0,ncol(M),ncol(oligo.data))

     # non-connected ones
     t1 <- which(rowSums(A>0)==1)
     H.nonneg[t1,] <- t(M[,t1]>0) %*% oligo.data

     # connected: divide into subsets
     t2 <- setdiff(1:nrow(A),t1)

     cmp.ndx <- conComp(A[t2,t2])

     for (i in 1:max(cmp.ndx)){
       i1 <- which(cmp.ndx == i)

       if (length(i1)>1)
         H.nonneg[t2[i1],] <- fcnnls(A[t2[i1],t2[i1]],B[t2[i1],], pseudo = TRUE)$x
       else {
         H.nonneg[t2[i1],] <- t(M[,t2[i1]]>0) %*% oligo.data
       }
     }
   } else {
     H.nonneg <- fcnnls(t(M) %*% M, t(M) %*% oligo.data, pseudo = TRUE)$x
   }

   colnames(H.nonneg) <- colnames(oligo.data)
   rownames(H.nonneg) <- colnames(M)

   return(H.nonneg)

}

