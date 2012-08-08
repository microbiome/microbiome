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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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


#' Description: level2-level1 mappings
#'
#' Arguments:
#'   @param l2 level2 phylotypes
#'   @param oligomap oligomap
#'
#' Returns:
#'   @return mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

level2TOlevel1 <- function (l2, oligomap) {

   # Check which L1/L2 id type used	       
   lev1 <- intersect(c("level 1", "level.1", "L1"), colnames(oligomap))
   lev2 <- intersect(c("level 2", "level.2", "L2"), colnames(oligomap))

   omap <- oligomap[match(as.character(l2), oligomap[[lev2]]), c(lev2, lev1)]
   omap[[lev2]] <- factor(omap[[lev2]])
   omap[[lev1]] <- factor(omap[[lev1]])

   omap
}

#' Description: L2-species mappings
#'
#' Arguments:
#'   @param l2 l2
#'   @param oligomap oligomap
#'
#' Returns:
#'   @return L2-species mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

L2.species <- function(l2, oligomap) {
  unique(as.character(oligomap$species[oligomap$level.2 == l2]))
}

#' Description: L2/L1-species mappings. For L1/L2 groups, list corresponding species.
#'
#' Arguments:
#'   @param level level
#'   @param oligomap oligomap
#'   @param groups return selected groups only
#'
#' Returns:
#'   @return level-species mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

level.species <- function(level, oligomap, groups = NULL) {
  res <- split(oligomap$species, oligomap[[level]])	      
  res <- lapply(res, unique)
  if (is.null(groups)) { groups <- names(res) }
  res[groups]
}




#' Description: species-levels mappings
#'
#' Arguments:
#'   @param spec species
#'   @param oligomap oligomap
#'
#' Returns:
#'   @return species-levels mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

species2levels <- function (spec, oligomap) {
   # Check which L1/L2 id type used	       
   lev1 <- intersect(c("level 1", "level.1", "L1"), colnames(oligomap))
   lev2 <- intersect(c("level 2", "level.2", "L2"), colnames(oligomap))
	       
   omap <- oligomap[match(as.character(spec), oligomap$species), c("species", lev2, lev1)]
   omap[["species"]] <- factor(omap[["species"]])
   omap[[lev1]] <- factor(omap[[lev1]])
   omap[[lev2]] <- factor(omap[[lev2]])

   omap
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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

#' sampleReproducibility function
#' @param con MySQL connection
#' @param sampleID sampleID
#'
#' @return TBA
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

sampleReproducibility <- function(con, sampleID) {
   if (phyloarrayConnection(con)) {
      if ((!(is.vector(sampleID) & length(sampleID)==1)) | (!is.character(sampleID))) {
         stop("Argument 'sampleID' must be a character vector of length 1")
      }
      extr <- fetch.extractions(con,condition=list(list(field='sampleID',value=sampleID)))
      extr <- extr[extr$isDiscarded==0 & extr$noReproCheck==0 & extr$noSampleNormalisation==0,]
      
   }
}
