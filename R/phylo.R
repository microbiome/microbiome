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


#' levelmap
#' 
#' Description: map phylotypes between hierarchy levels
#'
#' Arguments:
#'   @param phylotypes phylotypes to convert; if NULL then considering all phylotypes in the oligomap
#'   @param level.from level.from
#'   @param level.to level.to
#'   @param oligomap oligomap
#'
#' Returns:
#'   @return mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

levelmap <- function (phylotypes = NULL, level.from, level.to, oligomap) {

  oligomap <- polish.oligomap(oligomap)

  if (is.null(phylotypes)) {
    phylotypes <- as.character(unique(oligomap[[level.from]]))
  }

  if (level.from == "species" && level.to == "L2") {
    sl <- species2levels(phylotypes, oligomap)[c("species", "L2")]
  }	 

  if (level.from == "species" && level.to == "L1") {
    sl <- species2levels(phylotypes, oligomap)[c("species", "L1")]
  }	 

  if (level.from == "L2" && level.to == "L1") {
    sl <- level2TOlevel1(phylotypes, oligomap)[,2]
  }	 

  if (level.from == "L2" && level.to == "species") {
    sl <- list()
    for (pt in phylotypes) {
      sl[[pt]] <- as.character(unique(oligomap[oligomap[["L2"]] == pt, "species"]))
    }
  }	 

  if (level.from == "L1" && level.to == "L2") {
    sl <- list()
    for (pt in phylotypes) {
      sl[[pt]] <- as.character(unique(oligomap[oligomap[["L1"]] == pt, "L2"]))
    }
  }	 

  if (level.from == "L1" && level.to == "species") {
    sl <- list()
    for (pt in phylotypes) {
      sl[[pt]] <- as.character(unique(oligomap[oligomap[["L1"]] == pt, "species"]))
    }
  }	 

  sl

}



#' Phylogeneticenrichments
#'
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
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
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

#' retrieve.probesets
#' 
#' Description: List probes for each probeset
#'
#' Arguments:
#'   @param oligomap data.frame with oligo - phylotype mapping information
#'   @param level phylotype level for probesets
#' Returns:
#'   @return A list. Probes for each phylotype.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

retrieve.probesets <- function (oligomap, level = "species", name = NULL) {

  # If name not given, pick all
  if (is.null(name)) { name <- unique(as.character(oligomap[[level]])) }

  phylo <- oligomap[oligomap[[level]] %in% name,]

  if (is.factor(phylo[[level]])) {
    phylo[[level]] <- droplevels(phylo[[level]])
  }

  phylo.list <- split(phylo, phylo[[level]])
  probesets <- lapply(phylo.list, function(x) {as.character(unique(x$oligoID))})
  names(probesets) <- names(phylo.list)

  probesets

}
