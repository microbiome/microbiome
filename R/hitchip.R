#' Description: phylogenetic distance matrix for HITChip
#'
#' Arguments:
#'   @param level phylogeny level
#'   @param f.phylomap phylogenetic map with a column for each level, and one column for the ID which matches with the sequence similarity matrix in f. phylodist file
#'   @param f.phylodist probe sequence similarity matrix 
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{jarkko.salojarvi@@helsinki.fi}
#' @keywords utilities

hitchip.phylodistance <- function (level, f.phylomap, f.phylodist) {
  		     
  phylogeny <- read.csv(f.phylomap, sep = "\t", header = TRUE)
  tab <- read.csv(f.phylodist, sep = "\t", header = TRUE, row.names = 1)
  tab <- tab[colnames(tab), colnames(tab)]

  # Get all distances for each phylotype at given level
  pts <- unique(phylogeny[[gsub(" ", "", level)]])
  pd <- array(NA, dim = c(length(pts), length(pts)))
  rownames(pd) <- colnames(pd) <- pts
  for (pt1 in pts) {
    ids1 <- as.character(phylogeny$reference[which(phylogeny[[gsub(" ", "", level)]] == pt1)])
    for (pt2 in pts) {
      ids2 <- as.character(phylogeny$reference[which(phylogeny[[gsub(" ", "", level)]] == pt2)])
      # Use mean of pairwise phylogenies as the summary distance
      # for this level to the other levels
      if (length(ids1)>0 && length(ids2)>0) {
        mat <- as.matrix(tab[ids1, ids2], nrow = length(ids1))
        pd[pt1, pt2] <- 1 - mean(mat)
        pd[pt2, pt1] <- pd[pt1, pt2]
      }
    }
  }

  pd
}




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


#' Description: map phylotypes to upper hierarchy level
#'
#' Arguments:
#'   @param phylotypes phylotypes
#'   @param level.from level.from
#'   @param level.to level.to
#'   @param oligomap oligomap
#'
#' Returns:
#'   @return mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

levelmap <- function (phylotypes, level.from, level.to, oligomap) {

  levs3 <- c("species")	   
  levs2 <- c("level 2", "level.2", "L2")
  levs1 <- c("level 1", "level.1", "L1")

  lev3 <- intersect(levs3, colnames(oligomap))
  lev2 <- intersect(levs2, colnames(oligomap))
  lev1 <- intersect(levs1, colnames(oligomap))

  if (level.from %in% levs2) {level.from <- lev2}
  if (level.to %in% levs2) {level.to <- lev2}
  if (level.to %in% levs1) {level.to <- lev1}

  if (level.from == "species" && level.to == lev2) {
    sl <- species2levels(phylotypes, oligomap)[c("species", lev2)]
  }	 

  if (level.from == "species" && level.to == lev1) {
    sl <- species2levels(phylotypes, oligomap)[c("species", lev1)]
  }	 

  if (level.from == lev2 && level.to == lev1) {
    sl <- level2TOlevel1(phylotypes, oligomap)[,2]
  }	 

  if (level.from == lev2 && level.to == "species") {
    sl <- list()
    for (pt in phylotypes) {
      sl[[pt]] <- as.character(unique(oligomap[oligomap[[lev2]] == phylotypes, "species"]))
    }
  }	 

  if (level.from == lev1 && level.to == lev2) {
    sl <- oligomap[oligomap[[lev1]] == phylotypes, c(lev1, lev2)]
    sl <- sl[!duplicated(sl),2]
  }	 

  if (level.from == lev1 && level.to == lev3) {
    sl <- oligomap[oligomap[[lev1]] == phylotypes, c(lev1, lev3)]
    sl <- sl[!duplicated(sl),2]
  }	 

  sl

}



#' Description: simple background correction
#'
#' Arguments:
#'   @param dat data matrix
#'
#' Returns:
#'   @return background corrected data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

remove.background <- function (dat) {
  
  det.th <- estimate.min.threshold(dat)
  # Set observations below detection threshold to zero and shift other values to start from zero
  dat <- dat - det.th
  dat[dat<0] <- 0

  dat
} 


#' Description: determine threshold for bg correction
#'
#' Arguments:
#'   @param dat data matrix (in approximately normal scale ie. logged)
#'
#' Returns:
#'   @return threshold value
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

estimate.min.threshold <- function (dat) {
  #estimate min threshold 
  DD <- density(as.numeric(unlist(dat)))
  #find mode
  noise_mode <- DD$x[[which.max(DD$y)]] # DD$x[which(diff(DD$y)<0)[1]]
  #compute sd of noise
  noise_sd <- sd(dat[dat < noise_mode])
  #threshold
  low.thresh <- noise_mode + 6*noise_sd

  low.thresh
}


#' Description: Computes the dependence of the higher phylogenetic
#' level groups of the given list x, which can be a vector/data frame of either
#' oligos, species, or level2 groups. Computes only the enrichments
#' if onlyEnrich=T (neglects the unexpected disappearences of the
#' groups in the given list). Outputs only the enrichments below the
#' given p-value threshold p.th.
#' Requires oligomap=the first 4 columns from phylogenyprofile,
#' origlevel=the level from which the list is, and maplevel= the
#' level to which the mapping should be made and the enrichments of
#' which computed.
#' Outputs a four component list: p-values from Fisher's tests,
#' input tables for tests, the actual tests outputs, and the full
#' phylogenetic mapping information for the given list.
#' Changes: Version 1 computed two-tailed p-values also for single
#' occurrences on origlevel.
#' v2: two-tailed p-values, but not for single occurrences
#' v3: one-tailed p-values (enrichments), if wanted, outputting only
#' enrichments under the given p-value
#'
#' Arguments:
#'   @param x x
#'   @param oligomap oligomap
#'   @param origlevel origlevel
#'   @param maplevel maplevel 
#'   @param onlyEnriched onlyEnriched
#'   @param p.th p-value threshold
#'
#' Returns:
#'   @return enrichments
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

Phylogeneticenrichments <- function(x, oligomap, origlevel = colnames(oligomap)[3], maplevel = colnames(oligomap)[1], onlyEnriched = T, p.th = 0.05)
{
  
  ## Convert to character vector
  x <- as.character(x)
    
  ## collect origlevel groups for all x items
  origlevel.split <- split(oligomap,oligomap[[origlevel]])
  x.origlevelGroups <- origlevel.split[x]
  
  ## vector (inX and inXind) showing which items from the origlevel of
  ## oligomap are in the x list
  origlevel.ugroups <- names(origlevel.split)
  inXind <- match(x, origlevel.ugroups)
  inX <- rep(F,length(origlevel.ugroups))
  inX[inXind] <- T

  ## Collect the full phyloinfo
  phyloM <-  x.origlevelGroups[[1]]
  nolGroups <- length(x.origlevelGroups)
  if(nolGroups>1)
    for(i in 2:nolGroups){
      phyloM <- rbind(phyloM,x.origlevelGroups[[i]])
    }
  
  if(length(x)>2){
    
    ## compute enrichments of maplevel groups in x
    e <- c()
    pvals <- c()
    estimates <- c()
    tables <- list()
    maplevel.ugroups <- as.character(unique(phyloM[[maplevel]]))
    for(g in maplevel.ugroups){

      ## compute in which maplevel groups the given item occurs
      inMaplevelGroup <- unlist(lapply(origlevel.split,
                                         function(y){is.element(g,y[[maplevel]])
                                                   }))

      if(sum(inMaplevelGroup)>0){
        if(onlyEnriched)
          tmp <- try(fisher.test(inX, inMaplevelGroup, alternative="g"))
        else
          tmp <- try(fisher.test(inX, inMaplevelGroup, alternative="t"))
        if(tmp$p.value<p.th){
          e[g] <- list(tmp)
          pvals[g] <- tmp$p.value          
          tables[g] <- list(table(inX,inMaplevelGroup))
          estimates[g] <- tmp$estimate
        }
      }
    }

    ##Return enrichments i.e. Fisher's exact tests in table form

    if (length(pvals)>0)
       list(pvalues=t(rbind(pvals,estimates)), ##qvalues=t(t(qvalue(pvals))),
           tables=tables, tests=e, phyloMap=phyloM) 
    else
       list(pvalues=pvals, tables=tables, tests=e, phyloMap=phyloM) 

  } else
    list(pvalues=1, tables=NULL, tests=NULL, phyloMap=phyloM) 
}

