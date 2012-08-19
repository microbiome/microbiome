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

#' Description: Plot PCA with labels and explained variances
#'
#' Arguments:
#'   @param dPCA TBA
#'   @param main TBA
#'   @param xPC TBA
#'   @param yPC TBA
#'   @param showLabels TBA
#'   @param type TBA
#'   @param ... further parameters to be passed
#'
#' Returns:
#'   @return TBA
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plotPCA <- function(dPCA, main="PCA plot", xPC = 1, yPC = 2, showLabels = T, type = "n", ...){
   x11(width=11, height=11)
   dPCA.propVars <- round(dPCA$sdev^2/sum(dPCA$sdev^2)*100,digits=2)
   plot(dPCA$x[,xPC],dPCA$x[,yPC],type=type, main=main,
        xlab=paste("PC ",xPC," (",dPCA.propVars[xPC],"% of variance)",sep=""),
        ylab=paste("PC ",yPC," (",dPCA.propVars[yPC],"% of variance)",sep=""),...)
   if(showLabels)
     text(dPCA$x[,xPC],dPCA$x[,yPC]-0.02*max(dPCA$x[,yPC]),
          labels=rownames(dPCA$x), cex=0.7,...)
}



#' Description: 
#'
#' Arguments:
#'   @param dPCA TBA
#'   @param groups TBA
#'   @param main TBA
#'   @param xPC TBA
#'   @param yPC TBA
#'   @param type TBA
#'   @param legLoc TBA
#'   @param w TBA
#'   @param h TBA
#'   @param ... further parameters to be passed
#'
#' Returns:
#'   @return plot
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plotGroupPCA <- function(dPCA,groups,main="PCA plot",xPC=1, yPC=2,type="n",
                         legLoc="topright",w=11,h=11,...){
   x11(width=w, height=h)
   dPCA.propVars <- round(dPCA$sdev^2/sum(dPCA$sdev^2)*100,digits=2)
   plot(dPCA$x[,xPC],dPCA$x[,yPC],type=type, main=main,
        xlab=paste("PC ",xPC," (",dPCA.propVars[xPC],"% of variance)",sep=""),
        ylab=paste("PC ",yPC," (",dPCA.propVars[yPC],"% of variance)",sep=""),...)
   ##text(dPCA$x[,xPC],dPCA$x[,yPC]-0.02*max(dPCA$x[,yPC]),
   ##    labels=rownames(dPCA$x), cex=0.7,...)
   uGroups <- unique(groups)
   nG <- length(uGroups)
   gCols <- rainbow(nG)
   for(g in 1:nG){
     gFilter <- groups==uGroups[g]
     points(dPCA$x[gFilter,xPC],dPCA$x[gFilter,yPC],col=gCols[g], pch=19) 
   }
   legend(legLoc,legend=uGroups, fill=gCols)
 }

#' Description: Plot oligoprofile figures without data retrieval from the db
#'
#' Arguments:
#'   @param d data
#'   @param tax.level Taxonomic level
#'   @param metric metrics
#'   @param include.tree TBA
#'   @param fontsize font size
#'   @param figureratio figure ratio
#'   @param oligomap oligomap
#'   @param width width
#'   @param height height
#'   @param oligomapJN TBA
#'   @param plot.profile TBA
#'
#' Returns:
#'   @return plot
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plotOligoP <- function(d, tax.level='level1',
                        metric='correlation', include.tree=T,
                        fontsize=11, figureratio=15, oligomap=oligomapJN,
                        width=12, height=20, oligomapJN = NULL, 
			plot.profile = NULL){
 
   par(mar=c(0.1,0.1,0.1,0.1))
   x11(width=width, height=height)
   plot.profile(data=log(d), metric=metric,
                oligomap=oligomap, tax.level=tax.level,
                fontsize=fontsize, include.tree=include.tree,
                figureratio=figureratio)
}

#' Description: Plot Mean Profile 
#'
#' Arguments:
#'   @param data TBA
#'
#' Returns:
#'   @return barplot
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plotMeanProfile <- function(data){
  require(gplots)
  data.mean <- apply(data,1, mean)
  data.sd <- apply(data, 1, sd)
  barplot2(data.mean, horiz=T,
           plot.ci=T, ci.u=data.mean+data.sd, ci.l=data.mean-data.sd)
}


#' Description: Read data
#'
#' Arguments:
#'   @param ... No input by default
#'
#' Returns:
#'   @return Used for side effects
#'

#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read.in.data <- function(...) {
   tclvalue(rdir) <- paste(tclvalue(tkchooseDirectory(
                       title = "Read hitchip files from directory:")), "/",  sep="")
   tclvalue(wdir) <- paste(tclvalue(tkchooseDirectory(
                       title = "Save all the files produced by the script into directory:")), 
		         "/", sep = "")
   tclvalue(hitchip.info) <- tk_choose.files(
                               caption = "Select HITChip annotation file", 
			       multi = FALSE)
}

#' Generic method for submission of calculated data to the database 
#' @param x TBA
#'
#' @return TBA 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities


submitDB <- function(x) {
   if (!is.list(x)) {
      stop("Argument 'x' must be a list")
   }
   if (!is.null(x$dbcon)) {
      dbSendQuery(x$dbcon,'DROP TABLE IF EXISTS rinput')
   }
   else {
      stop("Missing 'dbcon' field in argument 'x'")
   }
   UseMethod("submitDB")
}

#' Delete previous outlier detection results for these samples. Also set
#' the reproducibility-check flag to 0 for all featureextractions belonging
#' to these samples.
#'
#' @param x TBA
#'
#' @return TBA 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

submitDB.sampleOutliers <- function( x ) {

   statement <- paste(
      "UPDATE featuremeasurement fm, featureextraction fe, hybridisation h ",
      "SET fm.isOutlier=0, fe.hasReproCheck=0 WHERE h.sampleID='", x$sampleID,"' ",
      "AND fm.extractionID=fe.extractionID ",
      "AND fe.hybridisationID=h.hybridisationID",
      sep='')

   rs <- dbSendQuery(x$dbcon, statement)

   # Put the outliers in a temporary table and transfer 
   # the data to the permanent tables
  if (nrow(x$spots) > 0) {
      dbWriteTable(x$dbcon,'Rinput',row.names=FALSE,x$spots,append=TRUE)
       dbSendQuery(x$dbcon,'ALTER TABLE rinput ADD INDEX (featureID)')
       dbSendQuery(x$dbcon,'ALTER TABLE rinput ADD INDEX (extractionID)')
       dbSendQuery(x$dbcon,'UPDATE featuremeasurement f, rinput r SET f.isOutlier=1 WHERE f.featureID=r.featureID and f.extractionID=r.extractionID')
  }
  
  NULL

}

#' panel.stability
#' @param x TBA
#' @param y TBA
#' @param scale TBA
#' @param ... parameters to pass
#'
#' @return TBA 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

panel.stability <- function (x, y, scale=c("linear","logarithmic"), ...) {
   scale <- match.arg(scale)
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   ## our x and y are already logarithms of signals
   if (scale=="logarithmic") {
      r <- abs(cor(x, y))
      m <- rsd(x, y)*100
   }
   if (scale=="linear") {
      r <- abs(cor(exp(x), exp(y)))
      m <- rsd(exp(x), exp(y))*100
   }
   txt <- paste('Pearson=',format(c(r, 0.123456789), digits=5)[1],sep='')
   txt <- paste(txt,paste('RSD=',format(c(m,0.123456789), digits=2)[1],'%',sep=''),sep="\n")
   text(0.5, 0.5, txt, cex=1.5)
}

#' Outlier detection based on replicate measurements of a sample
#'
#' This function takes all arrays on which a sample was analysed and
#' calculates, based on replicates of oligo's within and between arrays
#' which spots on which arrays are likely to be outliers.  The isOutlier
#' flag of the spot is set in the database if the spot is likely to be an
#' outlier.  The routine needs normalized data. If these are unavailable
#' for an array on which the sample was hybridised then the calculation
#' is skipped and a warning is issued in the R console.
#'
#' Note: sampleOutliers only considers probes in the class 'ssRNA', and
#' does not take into account hybridisations of which the 'isDiscarded'
#' flag is set, or featureExtractions of which the 'noReproCheck' is set.
#'
#'
#' @param con MySQL connection
#' @param sampleID sampleID
#' @param significance significance
#'
#' @return A sampleOutliers object which is a list containing a
#'           'spots' dataframe, a 'dbcon' connection to the phyloarray
#'           database from which the sample data were taken and a
#'           'sampleID' length 1 character vector containing the
#'           sampleID.
#' 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

# FIXME: to be removed?
sampleOutliers <- function (con, sampleID, significance=0.001) {
   if (phyloarrayConnection(con)) {
      if ((!(is.vector(sampleID) & length(sampleID)==1)) | (!is.character(sampleID))) {
         stop("Argument 'sampleID' must be a character vector of length 1")
      }
      out.extractionIDs <- c()
      out.featureIDs <- c()
      # Only ssRNA probes are selected, not control probes etc.
      statement <- paste(
         "SELECT o.oligoID, fm.extractionID, fm.featureID, log(spatNormSignal) AS signal ",
         "FROM hybridisation h, featureextraction fe, featuremeasurement fm, arrayfeature af, probe p, oligo o ",
         "WHERE h.sampleID='",sampleID,"' AND NOT h.isDiscarded ",
         "AND fe.hybridisationID=h.hybridisationID AND NOT fe.noReproCheck ",
         "AND fm.extractionID=fe.extractionID AND fm.featureID=af.featureID ",
         "AND af.probeID=p.probeID AND p.oligoID=o.oligoID AND o.class='ssRNA' ",
         "ORDER BY fm.extractionID, fm.featureID",
         sep='')
      rs <- dbSendQuery(con,statement)
      sig <- fetch(rs, n=-1)

      #sig <- quantnorm(sig) # in the original phyloarray
      sig <- ScaleProfile(sig, method = 'quant') # using R package function

      avgvar<-mean(tapply(sig$signal,sig$oligoID,var))
      haveoutlier <- tapply(sig$signal,sig$oligoID,outlierPvalue,avgvar)
      haveoutlier <- haveoutlier[haveoutlier<significance]
      if (length(haveoutlier)>0) {
         for (j in names(haveoutlier)) {
            #out.row <- sig[sig$oligoID==j,][minmax(sig[sig$oligoID==j,]$signal),]
            out.row <- sig[sig$oligoID==j,][ScaleProfile(sig[sig$oligoID==j,]$signal, method = 'minmax'),]
            out.extractionIDs <- c(out.extractionIDs,out.row$extractionID)
            out.featureIDs <- c(out.featureIDs,out.row$featureID)
         }
      }
      return(new("sampleOutliers",list(spots=data.frame(extractionID=out.extractionIDs,featureID=out.featureIDs), dbcon=con, sampleID=sampleID)))
   }
}


#' Description: Read HITChip file
#'
#' Arguments:
#'   @param hitchip.file hitchip file name
#'
#' Returns:
#'   @return matrix
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.hitchip2 <- function (hitchip.file) {

  amat <- read.csv(hitchip.file, sep = "\t", header = TRUE) 
  colnames(amat) <- gsub("\\.", "/", gsub("X", "", colnames(amat)))
  rownam <- amat[,1]
  amat <- amat[, -1]
  rownames(amat) <- rownam
  amat <- t(amat)

  # Remove spike-ins
  amat <- amat[, grep("Victivallis", colnames(amat), invert = T)]
  amat <- amat[, grep("Lentisphaerae", colnames(amat), invert = T)]

  amat

}


#' Description: read hitchip prof 008 v3.4
#'
#' Arguments:
#'   @param path TBA
#'   @param dataID TBA
#'
#' Returns:
#'   @return d
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.hitchip.prof008v3.4 <- function(path="../data/", dataID=""){

  ## Unweighted
  d <- list()
  d$o <- read.table(paste(path,dataID,"oligoprofile_008.tab", sep=""),
                    sep="\t",header=T, row.names=1)
  d$l1.sum <- read.table(paste(path,dataID,"level1_Sum_008.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.sum <- read.table(paste(path,dataID,"level2_Sum_008.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l1.log10.ave <- read.table(paste(path,dataID,"level1_log10Ave_008.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.log10.ave <- read.table(paste(path,dataID,"level2_log10Ave_008.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.log10.ave <- read.table(paste(path,dataID,"species_log10Ave_008.tab", sep=""),
                                  sep="\t",header=T, row.names=1)

  ## Weighted
  d$l1.sum.w <- read.table(paste(path,dataID,"level1_Sum_008_w.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.sum.w <- read.table(paste(path,dataID,"level2_Sum_008_w.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l1.log10.ave.w <- read.table(paste(path,dataID,"level1_log10Ave_008_w.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.log10.ave.w <- read.table(paste(path,dataID,"level2_log10Ave_008_w.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.log10.ave.w <- read.table(paste(path,dataID,"species_log10Ave_008_w.tab", sep=""),
                                  sep="\t",header=T, row.names=1)
  
  ## Specific
  d$l1.sum.s <- read.table(paste(path,dataID,"level1_Sum_008_s.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.sum.s <- read.table(paste(path,dataID,"level2_Sum_008_s.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l1.log10.ave.s <- read.table(paste(path,dataID,"level1_log10Ave_008_s.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.log10.ave.s <- read.table(paste(path,dataID,"species_log10Ave_008_s.tab", sep=""),
                                  sep="\t",header=T, row.names=1)

  return(d)
 }

#' Description: read hitchip prof 010
#'
#' Arguments:
#'   @param path TBA
#'   @param dataID TBA
#'
#' Returns:
#'   @return d
#'

#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @examples # usage example:
#' # d = read.hitchip.prof010("my.file.path")
#' # d = preproc.Hitchip(d) 
#' # mapply(doPCA,d,names(d))
#' # where doPCA is a function, for example
#' # doPCA=function(Data,filename){
#' #   pdf(paste(filename,"PCAplot.pdf",sep=""))
#' #   # screeplot:
#' #   plot(prcomp(t(Data)))
#' #   dev.off()
#' #}
#' #t.test example
#' # c.i a vector with class information
#' #t.test.res=apply(d$l2.log10.ave,1,function(x,y) t.test(x[y==1],x[y==2]),y=c.i)
#' #p.val.vec=sapply(t.test.res,function(x) x$p.value)
#' #estimate.vec=sapply(t.test.res,function(x) x$estimate)
#' #qvalue(p.val.vec)$qvalue
#' #which(p.val.vec<0.05)

read.hitchip.prof010 <- function(path="../data/", dataID = ""){

  ## Unweighted
  d <- list()
  d$o <- read.table(paste(path,dataID,"oligoprofile_010.tab", sep=""),
                    sep="\t",header=T, row.names=1)
  d$l1.sum <- read.table(paste(path,dataID,"level1_Sum_010.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.sum <- read.table(paste(path,dataID,"level2_Sum_010.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l1.log10.ave <- read.table(paste(path,dataID,"level1_log10Ave_010.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.log10.ave <- read.table(paste(path,dataID,"level2_log10Ave_010.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.log10.ave <- read.table(paste(path,dataID,"species_log10Ave_010.tab", sep=""),
                                  sep="\t",header=T, row.names=1)
  return(d)
}


#' Description: HITChip preprocessing function
#'
#' Arguments:
#'   @param x TBA
#'
#' Returns:
#'   @return x
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

preproc.Hitchip <- function (x) {
  x=lapply(x,function(z){ 
      if(length(grep("Victivallis",rownames(z)))>0){
         z=z[grep("Victivallis",rownames(z),invert=T),]} 
       return(z)})
  x=lapply(x,function(z){ 
      if(length(grep("Lentisphaerae",rownames(z)))>0){
         z=z[grep("Lentisphaerae",rownames(z),invert=T),]} 
       return(z)})
  x$o=x$o[!(rownames(x$o) %in% c("HIT 1503","HIT 1505","HIT 1506")),]
  #some oligos don't appear in phylogeny. Those can be removed by this:
  #x$o=x$o[rownames(x$o) %in% unique(PH.i[,5]),]
  #Here PH.i is Phylogenetic information phylogenyinfo_010.tab
  x$o=apply(x$o,2,log10)
  # log10 of relative abundances in sum data:
  x$l1.sum=apply(x$l1.sum,2,function(z) return(log10(z/sum(z))))
  # apply is a shorter way to write this: 
  # for (i in 1:ncol(x$l1.sum)){
  #   x$l1.sum[,i]=log10(x$l1.sum[,i]/sum(x$l1.sum[,i]))
  # }
  x$l2.sum=apply(x$l2.sum,2,function(z) return(log10(z/sum(z))))
  return(x)
}


#' Description: PCA dialogue window
#'
#' Arguments:
#'   @param d TBA
#'   @param d.i TBA
#'   @param write.dir Output directory
#'
#' Returns:
#'   @return tkpack
#'
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PCA.dialogue <- function(d, d.i, write.dir) {

  dlg <- tktoplevel()
  tkwm.deiconify(dlg)
  tkgrab.set(dlg)
  tkfocus(dlg)
  tkwm.title(dlg, "PCA")

  frm <- tkframe(dlg)
  frm.radio <- populate.radiobuttons(frm, title="HITChip Level", 
  	       				  var.names  = names(d), 
  	       				  var.values = names(d), 
					  var.init   = tclVar("l2.sum"))
  tkpack(frm,frm.radio$frame)
  label.set=c("None", colnames(d.i), "Full Annotation")
  frm.radio2=populate.radiobuttons(frm, title = "Labels",
				        var.names = label.set, 
					var.values = label.set,
					var.init = tclVar("None"))

  tkpack(frm,frm.radio2$frame,side="right")
  frm.check=populate.checkboxes(frm, title = "PC coordinates",
  				     var.names = c("PC1","PC2","PC3","PC4"),
				     var.init=c(1,1,0,0))

  tkpack(frm, frm.check$frame, side = "right")

  group.set=c("None",colnames(d.i)[2:ncol(d.i)])

  frm.check2 <- populate.radiobuttons(frm, title = "Grouping",
  	     				   var.names = group.set,
					   var.values = group.set,
					   var.init = tclVar("None"))

  tkpack(frm,frm.check2$frame, side = "right")

  button.run <- tkbutton(frm, text="Run",
                          command=function(){
     # gather all information from dialogue window
     t1 <- sapply(frm.check$var,tclvalue)
     t2 <- tclvalue(frm.check2$var)

     if (t2 == "None")
       group.by <- NULL
     else
       group.by <- d.i[[t2]]

     n1 <- tclvalue(frm.radio2$var)

     if (n1 == "None"){
       d.lab <- NULL
     } else {
       if (n1 == "Full Annotation")
         d.lab <- apply(d.i[,2:ncol(d.i)],1,paste,collapse="_")
       else
         d.lab <- as.character(d.i[[n1]])
    }

    # run PCA and plot
    runPCA(d,tclvalue(frm.radio$var),
		c1 = which(t1==1)[1], c2 = which(t1==1)[2],
		data.labels = d.lab,
		write.dir = write.dir, group.by = group.by)

    tkdestroy(dlg)
  })

  tkpack(frm, button.run, side = "right")

}


#' Description: low.level functions required by GUI
#'
#' Arguments:
#'   @param tt TBA
#'   @param title TBA
#'   @param var.names TBA 
#'   @param var.init TBA
#'
#' Returns:
#'   @return A list.
#'
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

populate.checkboxes <- function(tt, title, var.names, var.init){

 title.font <- tkfont.create(weight = "bold", size = 10)
 frm <- tkframe(tt, relief = "groove", borderwidth = 2)
 label.widget <- tklabel(frm, text = title, font = title.font)
 tkpack(label.widget, side = "top")
 cb <- list();
 var.val <- list()
 for (i in 1:length(var.names))
   var.val[[i]] <- tclVar(as.character(var.init[i]))
 names(var.val) <- var.names
 for (i in 1:length(var.val)){
   cb[[i]] <- tkcheckbutton(frm, text=names(var.val)[i],variable=var.val[[i]]);
   tkpack(cb[[i]], side="top", expand=1);
 }
 tkpack(frm, side = "top")
 
 return(list(frame = frm, var = var.val))

}




#' Description: Import Project From Phyloarray To R. This function imports a 
#' project's ln(data) from probeLevel from phyloarray mysql-database to R. 
#' Modified from avgSignalPerTaxon.sql Query 2 by JN
#' NOTE!!: Must be modified to average over the identical probes!! 
#' The desired projectName is given as an argument
#' as well as the algorithm versions.
#'
#' Arguments:
#'   @param projectName Project to import
#'   @param spatNormAlgv Spatial normalization algorithm version
#'   @param normAlgv Normalization algorithm version
#'   @param dbuser mysql database username
#'   @param dbpwd  mysql database password
#'   @param dbname mysql database name
#'
#' Returns:
#'   @return data
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

importProjectFromPhyloarrayToR <- function(projectName = 'ProvasI', 
			       spatNormAlgv = '1.1', normAlgv = '1.1', 
			       dbuser = NULL, dbpwd = NULL, dbname = NULL){

  # Intialize database connection
  drv <- dbDriver("MySQL")
  con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)

  # Get the projectID of the project
  projectID <- fetch(dbSendQuery(con,paste("SELECT p.projectID FROM project p WHERE p.projectName=\'",
                                         projectName,"\'",sep="")), n=-1)$projectID

  # Get all the sampleIDs from the given project
  statement <- paste(
                   "SELECT s.sampleID FROM sample s, hybridisation h, featureextraction f",
                   "WHERE s.sampleID=h.sampleID AND h.hybridisationID=f.hybridisationID",
                   "AND NOT h.isDiscarded AND NOT f.noReproCheck AND NOT f.noSampleNormalisation",
                   "AND s.projectID=\'",projectID,"\' AND spatNormAlgVersion>=",spatNormAlgv,
                   "GROUP BY sampleID HAVING COUNT(f.extractionID)>1")

  rs <- dbSendQuery(con,statement)
  fe <- fetch(rs, n = -1)
  print(length(fe))

  # Query and fetch quantile normalised data
  importRes <- c()
  if (length(fe)>0) {
    samples <- fe$sampleID
    for (i in samples) {
      message('Retrieving normalised sample "',i,'"\n')
      statement <- paste(
                         "SELECT probeName, p.oligoID, signal ",
                         "FROM probe p JOIN (",
                         " SELECT oligoID, avg(spatNormSignal) AS signal ",
                         " FROM project JOIN sample USING (projectID) ",
                         " JOIN hybridisation USING (sampleID) ",
                         " JOIN featureExtraction USING (hybridisationID) ",
                         " JOIN featureMeasurement USING (extractionID) ",
                         " JOIN arrayFeature USING (featureID) ",
                         " JOIN probe USING (probeID) ",
                         " WHERE projectName=\'",projectName,"\'",
                         " AND sampleID=\'",i,"\'",
                         " AND normalisationFinished ",
                         " AND NOT oligoID IS NULL ",
                         " AND NOT isDiscarded ",
                         " AND NOT noSampleNormalisation ",
                         " AND NOT isOutlier ",
                         " GROUP BY oligoID ",
                         ") AS os USING (oligoID) ",
                         "ORDER BY probeName; ", sep="")
                         
      rs <- dbSendQuery(con,statement)
      sig <- fetch(rs, n=-1)
      
      if (length(importRes)==0){
        importRes <- sig
        importRes$signal <- log(importRes$signal)
        names(importRes)[3] <- paste(as.character(i),sep="")
       } 
      else{
        importRes <- cbind(importRes,signal=log(sig$signal))
        names(importRes)[length(names(importRes))] <- paste(as.character(i),sep="")
      }
    }
  }

  dbDisconnect(con)

  return(importRes)
}




#' Description: plot scatter
#'
#' Arguments:
#'   @param sampleA TBA
#'   @param sampleB TBA
#'   @param d TBA
#'   @param cex TBA
#'
#' Returns:
#'   @return scatterplot 
#' 
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plotScatter <- function(sampleA, sampleB, d, cex=1){
  x11()
  plot(d[,sampleA],
       d[,sampleB], pch=".",
       xlab=sampleA, ylab=sampleB,
       main=paste("Between sample similarity, Pearson: ",
         round(cor(d[,sampleA],d[,sampleB]),2)))
  points(d[,sampleA],
         d[,sampleB], cex=cex, pch=".")
}


#' Description: TBA
#' 
#' Arguments:
#'   @param d.hitchip hitchip data matrix: phylogroups x samples
#'   @param annot annotations for samples, containing my.class, subjectID, sampleID, and perhaps time
#'   @param my.class specify the inspected class; 
#'   @param time if time is additionally used in modeling then annot needs to contain field named time
#'
#' Returns:
#'   @return TBA
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

hitchip.sig.lme <- function(d.hitchip, annot, my.class, time = NULL){
 
  # Define the class to be used in the analyses
  annot.tmp <- annot
  names(annot.tmp)[[which(names(annot.tmp) == my.class)]] <- "my.class"
  annot.tmp$my.class <- factor(annot.tmp$my.class)

  annot.tmp <- annot.tmp[which(annot.tmp$sampleID %in% colnames(d.hitchip)),]

  annot.tmp <- annot.tmp[apply(annot.tmp,1,function(z) sum(is.na(z)))==0,]

  coef.name <- "coefficients.1"

  d.hc <- d.hitchip[,annot.tmp$sampleID]

  if (length(levels(annot.tmp[[my.class]]))>1){
    if (length(levels(annot.tmp$time))>1){

      # Fit model for each phylogroup
      lme.list <- apply(d.hc, 1, function(y) { 
        annot.tmp$y <- y;
        M <- try( lme( y ~ time*my.class, random = ~1 | subjectID, data = annot.tmp), silent = T);

        if ( class(M)=="try-error" )
        M <- try( lme( y ~ time*my.class, random = ~1 | subjectID, data = annot.tmp, method = "ML"), silent = T);
         return(M)
        })

       retval <- lapply(lme.list, 
                   function(z){ 
                     if (class(z)!="try-error")
		       # FIXME: use package::glht to remove warnings in package build
                       p <- unlist(summary(glht(z,linfct=t(c(0,0,1,1))))$"test"[c("pvalues","coefficients")])
                     else
                       p <- NULL
                     return(p)})


        coef.name <- paste("+",names(coef(lme.list[[1]]))[c(3,4)],collapse="",sep="")

        retval <- lapply(retval,function(z){ names(z)[2]=coef.name; return(z)})

       } else { # No time information given

          lme.list <- apply(d.hc,1,function(y){annot.tmp$y=y; M=try(lm(y~my.class, data=annot.tmp),silent=T); return(M)})

          retval <- lapply(lme.list, function(z){ if (class(z)!="try-error")
	       # FIXME: use package::mcp to remove warnings in package build
               p=summary(glht(z,linfct=mcp(my.class ="Tukey")))$"test"[c("pvalues","coefficients")]
             else
               p=NA
             return(p)})
          coef.name="coefficients"
       }

    } else { # only a single class available
       lme.list <- apply(d.hc,1,function(y){annot.tmp$y=y; M=try(lme(y~time, random= ~1 | subjectID,data=annot.tmp),silent=T);
        if(class(M)=="try-error")
           M=try(lme(y~time, random= ~1 | subjectID, data=annot.tmp, method="ML"),silent=T);
        return(M)})
       retval=lapply(lme.list, function(z){ if (class(z)!="try-error")
            p=data.frame(pvalues=anova(z)["time",4], coefficients=coef(z)[1,2]) 
          else
            p=NULL
          return(p)})
          coef.name="coefficients"
    }

    # Models have been fitted, process:

    sig.res <- list(p.val=as.data.frame(lapply(retval,function(z) return(z["pvalues"]))),
     coef=as.data.frame(lapply(retval,function(z) return(z[coef.name]))))
    colnames(sig.res$coef)=rownames(d.hc)

    require(qvalue)
    if(min(dim(sig.res$p.val))==1){
       Q=try(qvalue(sig.res$p.val)$qvalue,silent=T)
       if(class(Q)=="try-error")
         Q=p.adjust(unlist(sig.res$p.val),method="BH")
       sig.res$q.val=Q
    } else { 
      sig.res$q.val=apply(sig.res$p.val,2,function(z){ Q=try(qvalue(z)$qvalue,silent=T)  

        if(class(Q)=="try-error"){
          Q <- p.adjust(z,method="BH")} 
          return(Q)
        })
    }
  
  return(sig.res)

}


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
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
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




#' Description: Load/install necessary packages and check OS
#'
#' Arguments:
#'
#' Returns:
#'   @return operating system string
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities


check.dependencies <- function () {

  ## determine if using windows or mac/linux, based on the path style
  if(strsplit(Sys.getenv()["R_HOME"],split="")[[1]][1]=="/"){
    os <- "unix"
  } else {
    os <- "win"
  }

  os
}



#' Description: simple background correction
#'
#' Arguments:
#'   @param dat data matrix
#'
#' Returns:
#'   @return background corrected data matrix
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

remove.background <- function (dat) {
  
  det.th <- estimate.min.threshold(dat)

  # Set observations below detection threshold to zero and shift other 
  # values to start from zero
  dat <- dat - det.th
  dat[dat<0] <- 0

  dat
} 


#' Description: list corresponding entities
#' 
#' Arguments:
#'   @param level level
#'   @param phylo phylo
#'
#' Returns:
#'   @return entities
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

list.entities <- function (level, phylo) {

  as.character(unique(phylo[[level]]))

}

#' Description: Returns the difference between the mean log10 signals of specific
#'              and unspecific probes per species
#'
#' Arguments:
#'   @param data TBA
#'   @param pi TBA
#'
#' Returns:
#'   @return TBA
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimateSpeciesPresence <- function(data, pi){
  
  pi.species <- split(pi, pi$species)

  est <- sapply(pi.species, function(x){
    spec <- x$nSpeciesPerOligo==1

    speciesd <- data[x$oligoID,]
    specid <- x$nSpeciesPerOligo==1
    uspecid <- x$nSpeciesPerOligo>1
    nSpecP <- sum(specid)
    nUspecP <- sum(uspecid)

    if(nSpecP>0 & nUspecP>0){
      statistic <- apply(speciesd,2,function(y){
        return(mean(log10(y[specid]), na.rm=T) - mean(log10(y[uspecid]), na.rm=T))
      })
    }else{
      if(nSpecP>0)
        statistic <- rep(5, dim(speciesd)[2])
      else
        statistic <- rep(-5, dim(speciesd)[2])
    }
    return(c(statistic, nSpecP = nSpecP, nUspecP = nUspecP))
  })

  return(est)
}


#' Description: read hitchip 
#'
#' Arguments:
#'   @param path Source path
#'   @param dataID data ID
#'
#' Returns:
#'   @return d
#'
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.hitchip <- function(path = "../data/", dataID = ""){
  d <- list()
  d$o <- read.table(paste(path,dataID,"oligoprofile.tab", sep=""),
                    sep="\t",header=T, row.names=1)

  d$l1.log10.sum <- read.table(paste(path,dataID,"level 1_log10_sum.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.log10.sum <- read.table(paste(path,dataID,"level 2_log10_sum.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.log10.sum <- read.table(paste(path,dataID,"species_log10_sum.tab", sep=""),
                                  sep="\t",header=T, row.names=1)

  d$l1.log10.ave <- read.table(paste(path,dataID,"level 1_log10_ave.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.log10.ave <- read.table(paste(path,dataID,"level 2_log10_ave.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.log10.ave <- read.table(paste(path,dataID,"species_log10_ave.tab", sep=""),
                                  sep="\t",header=T, row.names=1)

  d$l1.sum <- read.table(paste(path,dataID,"level 1_sum.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.sum <- read.table(paste(path,dataID,"level 2_sum.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.sum <- read.table(paste(path,dataID,"species_sum.tab", sep=""),
                                  sep="\t",header=T, row.names=1)

  d$l1.ave <- read.table(paste(path,dataID,"level 1_ave.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$l2.ave <- read.table(paste(path,dataID,"level 2_ave.tab", sep=""),
                             sep="\t",header=T, row.names=1)
  d$species.ave <- read.table(paste(path,dataID,"species_ave.tab", sep=""),
                                  sep="\t",header=T, row.names=1)

  try(d$o.divs <- read.table(paste(path,dataID,"oligodiversities.tab", sep=""),
                                  sep="\t",header=T, row.names=1))
  try(d$s.divs <- read.table(paste(path,dataID,"speciesdiversities.tab", sep=""),
                                  sep="\t",header=T, row.names=1))
  return(d)
  }


#' Description: get background parameters
#'
#' Arguments:
#'
#' Returns:
#'   @return TBA
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

get.bkg.params <- function(){
  tt <- tktoplevel()
  tkwm.title(tt, "Background Subtraction")

  frm <- populate.radiobuttons(tt,title="Background Subtraction Method",var.names=c("min. 500 oligos","2*sd bkg intensity","6*sd bkg intensity","none"),var.values=c("min. 500 oligos","2*sd bkg intensity","6*sd bkg intensity","none"),var.init=tclVar("2*sd bkg intensity"))

  frm2 <- tkframe(tt)

  frm2.up <- populate.radiobuttons(frm2,title="Handling of negatives in ave data",var.names=c("Keep negative values","Set negatives to zero"),var.values=c("TRUE","FALSE"),var.init=tclVar("FALSE"))

  frm2.down <- tkframe(frm2)
  button.OK <- tkbutton(frm2.down, text="OK", command=function(){
    tkdestroy(tt)
  })

  tkpack(frm2.down, button.OK)
  tkpack(frm2.up$frame,frm2.down,side="top",pady=5) 
  tkpack(frm2,frm$frame,side="left",padx=5) 
  tkpack(frm2)
  tkwait.window(tt) 
  return(list(method=tclvalue(frm$var),keep.neg=as.logical(tclvalue(frm2.up$var))))
}




#' 2D Spatial normalisation of microarrays
#' @param data dataframe with positional and signal information
#' @param x name of the x column in the dataframe
#' @param y name of the y column in the dataframe
#' @param signal name of the signal column in the dataframe
#' @param method method name, currently only "loess" is supported
#' @param span span
#' @param degree degree
#' @param family family
#' @param subset logical expression indicating rows to use for normalisation
#'
#' @return TBA
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

spatnorm <- function(data, x="x", y="y", signal="signal", method=c("loess"), span=0.03, degree=2, family="symmetric", subset=1==1) {
   method <- match.arg(method)
   # data integrity testing
   data <- as.data.frame(data)
   if (!is.data.frame(data)) {
      stop("Input should be a matrix or dataframe")
   }
   if (!(length(setdiff(c(x,y,signal),names(data)))==0)) {
    stop(paste("Missing columns in data:",paste(setdiff(c(x,y,signal),names(data)),collapse=" and ")))
   }
   # calculation
   fit <- loess(as.formula(paste(signal,"~",x,"*",y)),data,span=span,degree=degree,normalize=FALSE,family=family,control=loess.control(iterations=8,cell=0.07),subset=subset)
   pr.signal <- predict(fit,data)
   corr.signal <- data[[signal]] - pr.signal + min(pr.signal)
   corr.signal[corr.signal<0.1]=0.1
   outp <- list(fit=fit,pr.signal=pr.signal,corr.signal=corr.signal)
   invisible(outp)
}


#' Fetch measurements associated with one or more extractions from a phyloarray database
#' 
#' The function fetches records from the \emph{featureMeasurements}
#' table, including some fields from the array, hybridisation and
#' featureExtraction tables.
#' The \code{which.data} argument indicates which data to fetch from the
#' database, where \code{raw} will fetch the raw data (\code{FGsignal}),
#' \code{rawBGcor} will fetch the raw background corrected data
#' (\code{FGsignal-BGsignal}), \code{spatnorm} will fetch the spatially
#' normalized data and \code{qnorm} will fetch the sample-quantile
#' normalized data.
#' The \code{transformation} argument indicates which transformation to
#' apply to the data indicated by the \code{which.data} argument. Only
#' one type of transformation is allowed, where "\code{none}" will not
#' apply a transformation, "\code{avg}" will take the values averaged
#' over oligo's and featureExtractions including standard deviation,
#' "\code{log}" will take the 10-base logarithm of te data and
#' "\code{avglog}" will return the oligo and featureExtraction averaged
#' values of the logarithm of the data, including standard deviation.
#'
#' @param con MySQL connection
#' @param extrList An extractionList object obtained from fetch.extractions
#' @param which.data One or more of raw, rawBGcor, spatnorm and qnorm
#' @param transformation One of none, avg, log or avglog
#'
#' @return A dataframe with the selected measurements 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

fetch.measurements <- function (con, extrList, 
		   which.data = c("raw","rawBGcor","spatnorm","qnorm"), 
		   transformation = "none") {

   if (phyloarrayConnection(con)) {

      # Only one type of tranformation allowed
      tr <- match.arg(transformation,c("none","avg","log","avglog"),several.ok=FALSE)
      which.data <- match.arg(which.data, several.ok = TRUE)
      what <- list(none=list(id=c("af.featureID","oligoID","fm.isOutlier"),
                              raw="FGsignal AS raw",
                              rawBGcor="(FGsignal-BGsignal) AS rawBGcor",
                              spatnorm="spatNormSignal",
                              qnorm="normSignal"),

                   avg = list("id"="oligoID",
                              raw="avg(FGsignal) AS avgSignal, std(FGsignal) AS stdevSignal",
                             "rawBGcor"="avg(FGsignal-BGsignal) AS avgSignalBGcor, std(FGsignal-BGsignal) AS stdevSignalBGcor",
                             "spatnorm"="avg(spatNormSignal) AS avgSpatNormSignal, std(spatNormSignal) AS stdevSpatNormSignal",
                             "qnorm"="avg(normSignal) AS avgNormSignal, std(normSignal) AS stdevNormSignal"),

                  log = list("id"=c("af.featureID","oligoID","fm.isOutlier"),
                             "raw"="log(FGsignal) AS lgRaw",
                             "rawBGcor"="log(FGsignal-BGsignal) AS lgRawBGcor",
                             "spatnorm"="log(spatNormSignal) AS lgSpatNormSignal",
                             "qnorm"="log(normSignal) AS lgNormSignal"),

                  avglog = list("id"="oligoID",
                                "raw"="avg(log(FGsignal)) AS avgLgSignal, std(log(FGsignal)) AS stdevLgSignal",
                                "rawBGcor"="avg(log(FGsignal-BGsignal)) AS avgLgSignalBGcor, std(log(FGsignal-BGsignal)) AS stdevLgSignalBGcor",
                                "spatnorm"="avg(log(spatNormSignal)) AS avgLgSpatNormSignal, std(log(spatNormSignal)) AS stdevLgSpatNormSignal",
                                "qnorm"="avg(log(normSignal)) AS avgLgNormSignal, std(log(normSignal)) AS stdevLgNormSignal")          
                  )

      fields <- paste(c(c("sampleID","fe.extractionID"),
                       what[[tr]]$id,unlist(what[[tr]][which.data])), 
		       collapse=", ")

      grouping <- ifelse((tr=="avg" | tr=="avglog"),"GROUP BY oligoID, extractionID","")

      condition <- expandCondition(list(list(field='extractionID',value=extrList$extractions$extractionID)))

      stm <- paste("SELECT",fields,"FROM hybridisation h ",
                   "JOIN featureExtraction fe USING (hybridisationID)",
                   "JOIN featureMeasurement fm USING (extractionID)",
                   "JOIN arrayFeature af USING (featureID)",
                   "JOIN probe p USING (probeID)",condition, grouping)

      rs <- dbSendQuery(con, stm)
      data <- fetch(rs, n=-1)
      return(data) 

   }   
}


#' Relative standard deviation using paired measurements
#'
#'  Calculates an unbiased estimate of the the relative squared
#'  difference of x.  Although most RSD calculations use the following
#'  formula: mean(abs(x_i-y_i)/(x_i+y_i)) we here use:
#'  mean(sqrt(pi)*abs(x_i-y_i)/(x_i+y_i)) The latter is an unbiased
#'  estimate of the relative standard deviation (sigma/mu) when data are
#'  normal distributed.
#'
#' @param x dataframe, matrix or vector
#' @param y dataframe, matrix or vector compatible with x
#'
#' @return Returns a matrix of rsd values. 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

rsd <- function(x, y = NULL) {
   if (is.data.frame(y)) {
      y <- as.matrix(y)
   } else {
      stopifnot(is.atomic(y))
   }
   if (is.data.frame(x)) {
      x <- as.matrix(x)
   } else {
      stopifnot(is.atomic(x))
      if (!is.matrix(x)) {
         if (is.null(y))
            stop("Supply both 'x' and 'y' or a matrix-like 'x'")
         x <- as.vector(x)
      }
   }
   if (!is.null(y)) {
      x <- cbind(x,y)
   }
   RSD <- matrix(,nrow=dim(x)[2],ncol=dim(x)[2],dimnames=list(colnames(x),colnames(x)))
   for (i in 1:dim(x)[2]) {
      for (j in 1:dim(x)[2]) {
         RSD[i,j]=mean(sqrt(pi)*abs(x[,i]-x[,j])/abs(x[,i]+x[,j]),na.rm=TRUE)
      }
   }
   #if (dim(x)[2]==2) {
   #   RSD <- RSD[1,2]
   #}
   return(RSD)
}



#' Statistical test for outliers in a vector
#'
#' Function for determining the significances (p-values) of outliers in vector.
#' Uses chisq.out.test from the package 'outliers'.
#'
#' Arguments:
#'  @param x numerical vector
#'  @param avgvar the estimated variance of the distribution
#' Returns:
#'  @return Returns the p-value of the outlier test 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

outlierPvalue <- function (x, avgvar) {
   test <- chisq.out.test(x, variance = avgvar)
   return(test$p.value)
}


