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

# ------------------------------------------------------------


#' Description: Cross-correlate input variables
#'              
#' Arguments:
#'   @param annot annotation matrix: samples x features
#'   @param dat data matrix: samples x features
#'   @param method association method (pearson, spearman for continuous; categorical for discrete
#'   @param qth q-value threshold for included features 
#'   @param cth correlation threshold
#'   @param order order the results
#'   @param n.signif mininum number of significant correlations for each element
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with cor, pval, qval, N
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

cross.correlate <- function(annot, dat, method = "pearson", qth = NULL, cth = NULL, order = FALSE, n.signif = 0, verbose = TRUE, mode = "list"){

  x <- annot # numeric or discrete
  y <- dat # numeric

  if (is.null(colnames(y))) {colnames(y) <- paste("column-", 1:ncol(y), sep = "")}

  xnames <- colnames(x)
  ynames <- colnames(y)
  qv <- nmat <- NULL

  # Rows paired.

  Pc <- matrix(NA, ncol(x), ncol(y))
  Cc <- matrix(NA, ncol(x), ncol(y))
  rownames(Cc) <- colnames(x)
  colnames(Cc) <- colnames(y)
  rownames(Pc) <- colnames(x)
  colnames(Pc) <- colnames(y)


  # Ensure the correct formats

  if (ncol(x) == 1) {
    x <- matrix(t(apply(x, 1, as.numeric)), nrow = length(x))
  } else {
    if (!method == "categorical") {
      x <- matrix(t(apply(x, 1, as.numeric)), nrow(x))
    } else {
      x <- matrix(t(apply(x, 1, as.factor)), nrow(x))
    }
  }
  colnames(x) <- xnames

  if (ncol(y) == 1) {
    y <- matrix(t(apply(y, 1, as.numeric)), nrow = length(y))
  } else {
      y <- matrix(t(apply(y, 1, as.numeric)), nrow(y))
  }
  colnames(y) <- ynames
  
  # ----------------------------------------------------------------

  if (method %in% c("pearson","spearman")) {

    for (j in 1:ncol(y)){
      if (verbose) {print(j/ncol(y))}
      jc <- apply(x, 2, function (xi) {res <- cor.test(xi, y[, j], method = method, use = "pairwise.complete.obs"); c(res$estimate, res$p.value)})

      Pc[,j] <- jc[2,]        
      Cc[,j] <- jc[1,]        

    } 

  } else if (method == "bicor") {
    require(WGCNA)
    t1 <- bicorAndPvalue(x, y, use = "pairwise.complete.obs")
    Pc <- t1$p
    Cc <- t1$bicor
   } else if (method == "categorical") {  

      Cc <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
      rownames(Cc) <- colnames(x)
      colnames(Cc) <- colnames(y)
      nmat <- Pc <- Cc

      for (varname in colnames(x)) {

        for (lev in colnames(y)) {

          xvec <- x[, varname]
	  yvec <- y[, lev]
	  keep <- rowSums(is.na(cbind(xvec,yvec))) == 0	       
	  xvec <- xvec[keep]
	  yvec <- yvec[keep]
	
	  # Number of data-annotation samples for calculating the correlations
    	  n <- sum(keep)

    	  Cc[varname, lev] <- GKtau(xvec, yvec)

        }
      }
   } else if (method == "mi") {

      Cc <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
      rownames(Cc) <- colnames(x)
      colnames(Cc) <- colnames(y)
      nmat <- Pc <- Cc

      for (i in 1:ncol(x)) {
        for (j in 1:ncol(y)) {

          Cc[i,j] <- build.mim(cbind(x[,i], y[,j]), estimator = "spearman")[1, 2]
        }
      }
   }

   if (!all(is.na(Pc))) {

     # Corrected p-values
     require(qvalue)
     if (prod(dim(Pc)) >= 100) {
       qv <- qvalue(Pc)
       if (length(qv) == 1) { 
         qv <- NULL 
       } else {
         qv <- qv$qvalue
         rownames(qv) <- xnames
         colnames(qv) <- ynames
       }
     } else {
       cat("Too few p-values available, q-value calculation skipped-")
       qv <- NULL
     }

     rownames(Pc) <- xnames
     colnames(Pc) <- ynames

     rownames(Cc) <- xnames
     colnames(Cc) <- ynames
  }

   # Filter
   if (!is.null(qth) || !is.null(cth)) {

     # Filter by qvalues and correlations
     inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL

     if (!is.null(qth)) {
       inds1.q <- apply(abs(qv), 1, function(x) {sum(x < qth) >= n.signif}) 
       inds2.q <- apply(abs(qv), 2, function(x) {sum(x < qth) >= n.signif}) 
     }

     if (!is.null(cth)) {
       inds1.c <- apply(abs(Cc), 1, function(x) {sum(x > cth) >= n.signif})
       inds2.c <- apply(abs(Cc), 2, function(x) {sum(x > cth) >= n.signif})
     }

     if (!is.null(qth) && !is.null(cth)) {
       inds1 <- inds1.q & inds1.c
       inds2 <- inds2.q & inds2.c
     } else if (is.null(qth) && !is.null(cth)) {
       inds1 <- inds1.c
       inds2 <- inds2.c
     } else if (!is.null(qth) && is.null(cth)) {
       inds1 <- inds1.q
       inds2 <- inds2.q
     }

     Cmat <- as.matrix(0)
     # TODO: add also correlation filter, not only significance
     # Require each has at least n.signif. correlations
     if (sum(inds1) >= n.signif && sum(inds2) >=n.signif) {

       rnams <- rownames(Cc)[inds1]
       cnams <- colnames(Cc)[inds2]

       Cc <- matrix(Cc[inds1,inds2, drop = FALSE], nrow = sum(inds1))
       Pc <- matrix(Pc[inds1,inds2, drop = FALSE], nrow = sum(inds1))
       qv <- matrix(qv[inds1,inds2, drop = FALSE], nrow = sum(inds1))

       rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
       colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams

       if ( order && sum(inds1) >=2 && sum(inds2)>=2 ) { # Order in visually appealing order

         tmp <- Cc
         rownames(tmp) <- NULL
         colnames(tmp) <- NULL
         h <- heatmap(tmp, xlab = NULL, ylab = NULL, xaxt = 'n', yaxt = 'n'); dev.off()

         rnams <- rownames(Cc)[h$rowInd]
         cnams <- colnames(Cc)[h$colInd]
         Cc <- Cc[h$rowInd,h$colInd]
         Pc <- Pc[h$rowInd,h$colInd]
         qv <- qv[h$rowInd,h$colInd]

         rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
         colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams

       }

     } else {
       cat("No significant correlations with the given criteria\n")
       Cc <- Pc <- qv <- NULL
    }
   }

   res <- list(cor = Cc, pval = Pc, qval = qv, N = nmat)

   if (mode == "list") {
     return(res)     
   } else if (mode == "table") {
     ctab <- cbind(melt(res$cor), melt(res$qval)$value)
     colnames(ctab) <- c("X1", "X2", "correlation", "qvalue")
     o <- order(ctab$qvalue)
     ctab <- ctab[o, ]
     return(ctab)
   }
}



