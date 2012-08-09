#' levelmap
#' 
#' Description: map phylotypes between hierarchy levels
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

#' list.probes
#' 
#' FIXME: merge with retrieve probesets
#' 
#' Description: list probes corresponding to a given phylogenetic level
#' 
#' Arguments:
#'   @param name name
#'   @param level level
#'   @param phylo phylo
#'
#' Returns:
#'   @return probes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.probes <- function (name, level, phylo) {

  as.character(unique(phylo[which(phylo[[level]] == name), "oligoID"]))
  
}

#' retrieve.probesets
#' 
#' Description: List probes for each probeset
#'
#' Arguments:
#'   @param phylo data.frame with oligo - phylotype mapping information
#'   @param level phylotype level for probesets
#' Returns:
#'   @return A list. Probes for each phylotype.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

retrieve.probesets <- function (phylo, level = "species") {

  # phylo <- pruned16S
  phylo.list <- split(phylo, phylo[[level]])
  probesets <- lapply(phylo.list, function(x) {unique(x$oligoID)})	     
  names(probesets) <- names(phylo.list)

  probesets

}
