# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <leo.lahti@iki.fi>. All rights reserved.

# This file is a part of the microbiome R package

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' Description: PCA function
#'
#' Arguments:
#'   @param d TBA
#'   @param level TBA 
#'   @param c1 TBA
#'   @param c2 TBA
#'   @param data.labels Optional
#'   @param write.dir Output directory path
#'   @param group.by Optional
#'   @param ... Other parameters to be passed to the function
#'
#' Returns:
#'   @return A list.
#'
#' @export
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

runPCA <- function(d, level, c1, c2, data.labels = NULL, write.dir, group.by = NULL, ...){

  y <- as.matrix(d[[level]])
  mm <- prcomp(t(y),retx=T)
  mm.propVars <- round(mm$sdev^2/sum(mm$sdev^2)*100,digits=2)

  #screeplot
  plot(mm,main="Eigenvalues")

  pdf(paste(write.dir,"PCA_scree.pdf",sep=""))
  plot(mm,main="Eigenvalues")
  dev.off()
  if (is.null(group.by))
     col.vec="black"
  else{
     col.vec=as.numeric(as.factor(group.by))
     col.vec=colorRampPalette(c("Red", "Green","Blue"),space="rgb")(max(col.vec))[col.vec]
  }
  dev.new()
  
  plot(mm$x[,c1],mm$x[,c2],type="n", main="PCA",
       xlab=paste("PC ",c1," (",mm.propVars[c1],"% of variance)",sep=""),
       ylab=paste("PC ",c2," (",mm.propVars[c2],"% of variance)",sep=""),...)

  if (is.null(data.labels))
    points(mm$x[,c1],mm$x[,c2], cex=0.7,pch=20,col=col.vec)
  else
    text(mm$x[,c1],mm$x[,c2]-0.02*max(mm$x[,c2]),col=col.vec,labels=data.labels, cex=0.7,...)

  pdf(paste(write.dir,"PCA_12",sub(".","",level,fixed=T),".pdf",sep=""))
  plot(mm$x[,c1],mm$x[,c2],type="n", main="PCA",
       xlab=paste("PC ",c1," (",mm.propVars[c1],"% of variance)",sep=""),
       ylab=paste("PC ",c2," (",mm.propVars[c2],"% of variance)",sep=""),...)
  if (is.null(data.labels))
    points(mm$x[,c1],mm$x[,c2], cex=0.7,pch=20,col=col.vec)
  else
    text(mm$x[,c1],mm$x[,c2]-0.02*max(mm$x[,c2]),col=col.vec,labels=data.labels, cex=0.7,...)
  dev.off()

}




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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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



#' Description:  Computes the given lm model for each variable (row) in the 
#'               given matrix/data frame. Designed for HITchip matrices, 
#' 		 but is applicable to any other matrix also.
#'  		 NOTE: does not take log! 
#'
#' FIXME: merge with lm.matrix2
#'
#' Arguments:
#'   @param formula formula
#'   @param d.matrix data matrix
#'   @param d.info information 
#'   @param type information to return
#'
#' Returns:
#'   @return matrix
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities


lm.matrix <- function(formula, d.matrix, d.info, type = "coef") {

  require(limma)

  lm.res <- NULL

  for (v in colnames(d.matrix)) {

      df <- data.frame(list(d.info, y = as.vector(d.matrix[,v])))
      lmfit <- lm(formula, df)
      if (type == "coef") {
        lm.res <- rbind(lm.res, lmfit$coefficients)
      } else if (type == "qval") {
        lm.res <- c(lm.res, anova(lmfit)[["Pr(>F)"]][[1]])
      }
  }

  if (type == "coef") {
    rownames(lm.res) <- colnames(d.matrix)
  } else if (type == "qval") {
    names(lm.res) <- colnames(d.matrix)
    require(qvalue)
    lm.res[is.na(lm.res)] <- 1 # pvalue = 1 for missing vals
    lm.res <- qvalue(lm.res, pi0.method = "bootstrap")$qvalue
  }

  lm.res

}


#' Description:  Computes the given lm model for each variable (row) in the given matrix/data frame
#'  		 Designed for HITchip matrices, but is applicable to any other matrix also.
#'  		 NOTE: does not take log automatically!
#'
#' Arguments:
#'   @param formula formula
#'   @param d.matrix data matrix
#'   @param d.info information 
#'   @param contr TBA
#'   @param Cvar TBA
#'   @param fit.contrast 
#'
#' Returns:
#'   @return lmMatrixRes
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

lm.matrix2 <- function(formula, d.matrix, d.info, contr = NULL, Cvar = NULL, fit.contrast = NULL){ # na.action = na.fail)
  
  require(limma)

  vars <- rownames(d.matrix)
  lmMatrixRes <- list()

  cat("Computing lm for the following variables and data:\n")
  print(vars)
  str(d.info)
  str(t(d.matrix))
  
  for(v in vars){
    data <- d.info
    data$y <- as.vector(t(d.matrix[v,]))

    lmMatrixRes$lm[[v]] <- lm(formula, data)#, na.action)
    lmMatrixRes$anova[[v]] <- try(anova(lmMatrixRes$lm[[v]]))
    lmMatrixRes$summary[[v]] <- try(summary(lmMatrixRes$lm[[v]]))
    if (!is.null(Cvar)){
      lmMatrixRes$contrast[[v]] <- fit.contrast(lm(formula, data), Cvar, contr)
    }
  }
  lmMatrixRes
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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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



#' Description: Correlation distance. Computes the correlation 
#' distance between the _columns_ of the given matrix
#'
#' Arguments:
#'   @param x data matrix
#'
#' Returns:
#'   @return distance object
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

corDist <- function(x){
  as.dist((1-cor(x, use = "pairwise.complete")))
}


#' Description: wardClust is for heatmap function for Ward-clustering
#'
#' Arguments:
#'   @param x data matrix
#'
#' Returns:
#'   @return hclust object
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

wardClust <- function(x){hclust(x,method = "ward")}

#' Description: heatmap function for complete linkage clustering
#'
#' Arguments:
#'   @param x data matrix
#'
#' Returns:
#'   @return hclust object
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

completeClust <- function (x) {
  hclust(x, method = "complete")
}

#' Description: Perform pairwise comparison between factor levels
#'              
#' Arguments:
#'   @param x data matrix: samples x features
#'   @param y factor levels for the samples
#'   @param qth q-value threshold for included features 
#'   @param resdir Optional. If given, output files will be written in the result directory
#'
#' Returns:
#'   @return List with top findings from pairwise comparisons and their q-values
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

pairwise.comparisons <- function (x, y, qth = 0.05, resdir = NULL) {

  y <- factor(y)

  library(limma)

  # Remove categories with NA annotations and very small sample size
  keep <- !is.na(y) & (y %in% names(which(table(y) > 2)))
  y <- y[keep]
  x <- x[keep, ]

  top.findings <- list()
  top.findings.qvals <- list()

  # Compare each group to other groups
  for (i in 1:(length(unique(y))-1)) {

    qval.list <- list()

    for (j in (i+1):(length(unique(y)))) {
   
      ai <- as.character(unique(y)[[i]])
      aj <- as.character(unique(y)[[j]])

      inds1 <- which(y == ai)
      inds2 <- which(y == aj)

      lab <- factor(y[c(inds1, inds2)])
      ann <- data.frame(list(varname = lab))
      dat <- x[c(inds1, inds2),]

      # Here y refers to internal lm.matrix response variable(!)
      qvals <- lm.matrix(y ~ varname,  dat, ann, type = "qval")
      coefs <- lm.matrix(y ~ varname - 1,  dat, ann, type = "coef") # remove intercept

      nams <- names(sort(qvals))
      
      fc <- colMeans(x[inds1, nams]) - colMeans(x[inds2, nams])
    
      tab <- cbind(nams, coefs[nams, paste("varname", ai, sep = "")], coefs[nams, paste("varname", aj, sep = "")])
      colnames(tab) <- c("phylotype", ai, aj)
      tab <- as.data.frame(tab)
      tab$fold.change <- fc[nams]
      tab$qvalue <- qvals[nams]

      top.findings[[paste(ai, aj, sep = "-")]] <- subset(tab, qvalue < qth)
      qval.list[[aj]] <- qvals
      top.findings.qvals[[ai]] <- qval.list

      if (!is.null(resdir)) {

        resfile <- paste(resdir,"/limma-", ai, "-", aj, "-results.tab", sep = "")
        message(paste("Writing comparison statistics in: ",  resfile))
        write.table(tab, file = resfile, sep = "\t", quote = FALSE, row.names = FALSE)

      if (nrow(tab) > 0) {
    
        df <- as.data.frame(t(rbind(group = as.character(lab), t(dat)[rownames(tab)[1:min(nrow(tab), 9)],])))
        df$group <- factor(df$group)
        dfm <- reshape::melt.data.frame(df, id.vars = c("group"))
        dfm$value <- as.numeric(as.character(dfm$value))
        dfm$title <- apply(cbind(as.character(dfm$variable), round(as.numeric(as.character(tab$qvalue[match(dfm$variable, names(tab$qvalue))])),3)),1,function(x){paste(x, collapse = "/ qval: ")})
        p <- ggplot(dfm, aes(group, value)) + geom_boxplot() + facet_wrap(~title) + opts(title = "Top hits")
        pdfname <- paste(resdir,"/limma-", ai, "-", aj, "-tophit.boxplots.pdf", sep = "")
        print(pdfname)
        pdf(pdfname)
        print(p)
        dev.off()
  
      }
     }
    }
  }

  list(qval = top.findings.qvals, top = top.findings)

}

############################

check.wilcoxon <- function () {

  require(svDialogs)

  ## Open your tab file, Level 1&2 Sum_BGsub_Rel.contribution

  file <- choose.files(multi=F)
  oligo <- read.table(file,sep="\t", header=T, row.names=1)
  samples <- colnames(oligo)

  ## To select samples you can do that in 2 ways: select G1 and those samples not in G1 are G2, for which you would use:
  ## G2 <- samples[!(samples %in% G1)]
  ## or select G1 and select G2 (this is useful when you have multiple groups) (standard below)

  G1 <- tk_select.list(samples, multiple=T, title="Select samples for 1st group")
  G2 <- tk_select.list(samples, multiple=T, title="Select samples for 2nd group")

  levels <- rownames(oligo)

  M <- matrix(data=NA,length(levels),1)
  rownames(M) <- levels

  for (i in 1:length(levels)) {
	
	lvl <- levels[i]
	l.g1 <- oligo[lvl,G1]
	l.g2 <- oligo[lvl,G2]
	
	p<-wilcox.test(as.numeric(l.g1),as.numeric(l.g2))$p.value

	

	cat(lvl," p-value: ",p, "\n")

	M[i,1]<-p
	}

  M

  ## To Adjust P-values for Multiple Comparisons with Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
  cor.p <- p.adjust(M,method="BH") # Other methods can be applied here
  cor.p
}
