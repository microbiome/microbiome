# Copyright (C) 2006-2012 Douwe Molenaar, Janne Nikkilä, Leo Lahti, and Jarkko Salojarvi 
#
# Contact: <leo.lahti@@iki.fi>. All rights reserved.
#
# This file is a part of the microbiome R package
#
# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file has been modified from the previous phyloarray script 
# originally written by DM. 
#
# FIXME: Some of these function may be redundant or unnecessary, consider
# removal.


#' Description: Select projects to analyze
#' 
#' Arguments:
#'   @param con valid MySQL connection
#'   @param multi enable selection of multiple options
#'   @param condition TBA
#' Returns:
#'   @return vector of project names 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

choose.projects <- function (con, multi = TRUE, condition = NULL) {
   prjs <- fetch.projects(con, condition = condition)
   projects <- select.list(sort(prjs$projectName), multiple = multi, title = "Select studies:")
   prjs <- fetch.projects(con, condition = list(list(field = 'projectName', value = projects)))
   return(prjs)
}

######################################


#' Fetch projects from the phyloarray database
#' @param con MySQL connection
#' @param condition list of lists with field-value pairs
#'
#'  The function fetches complete records from the \emph{projects}
#' table. If \code{condition=NULL} it will fetch all records and if
#' \code{condition} is defined it will fetch only those records that
#' comply with the condition.  Condition must be a list of lists, each
#' having at least the fields \emph{field} and \emph{value}.  The
#' \emph{field} field must be a character vector of length 1 and
#' \emph{value} must be a vector.  Each of the values will be evaluated
#' as a optional value for the \emph{field} field.  This is an example:
#' "\code{list(list(field='projectName',value=c(A,B)))}" that will be
#' evaluated to the SQL condition "\code{WHERE (projectName='A' OR
#' projectName='B')}".  
#' @return A dataframe with the selected records from the projects table
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

fetch.projects <- function (con, condition = NULL) {
   if (phyloarrayConnection(con)) {
      stm <- paste("SELECT * FROM project", expandCondition(condition), sep='')
      rs <- dbSendQuery(con, stm)
      prjs <- fetch(rs, n=-1)
      return(prjs)
   }
}


#' Fetch samples from a phyloarray database
#'
#' The function fetches complete records from the \emph{samples}
#' table. If \code{condition=NULL} it will fetch all records and if
#' \code{condition} is defined it will fetch only those records that
#' comply with the condition.  Condition must be a list of lists, each
#' having at least the fields \emph{field} and \emph{value}.  The
#' \emph{field} field must be a character vector of length 1 and
#' \emph{value} must be a vector.  Each of the values will be evaluated
#' as a optional value for the \emph{field} field.  This is an example:
#' "\code{list(list(field='projectName',value=c(A,B)))}" that will be
#' evaluated to the SQL condition "\code{WHERE (projectName='A' OR
#' projectName='B')}".
#'
#' @param con MySQL connection
#' @param condition list of lists with field-value pairs
#'
#' @return A dataframe with the selected records from the samples table
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

fetch.samples <- function (con, condition = NULL) {
   if (phyloarrayConnection(con)) {
      stm <- paste("SELECT * FROM sample", expandCondition(condition), sep='')
      rs <- dbSendQuery(con, stm)
      smps <- fetch(rs, n=-1)
      return(smps)
   }
}


#' fetch.extractions
#' Fetch extractions from the phyloarray database.
#'
#' The function fetches complete records from the \emph{extractions}
#' table, including some fields from the array and hybridisation
#' tables. If \code{condition=NULL} it will fetch all records and if
#' \code{condition} is defined it will fetch only those records that
#' comply with the condition.  Condition must be a list of lists, each
#' having at least the fields \emph{field} and \emph{value}.  The
#' \emph{field} field must be a character vector of length 1 and
#' \emph{value} must be a vector.  Each of the values will be evaluated
#' as a optional value for the \emph{field} field.  This is an example:
#' "\code{list(list(field='projectName',value=c(A,B)))}" that will be
#' evaluated to the SQL condition "\code{WHERE (projectName='A' OR
#' projectName='B')}".
#'
#' @param con MySQL connection
#' @param condition list of lists with field-value pairs
#' @return An \code{extractionList} object consisting of a list with one entry called "extractions" which is a dataframe consisting of the records obtained from the database. 
#'
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities


fetch.extractions <- function (con, condition = NULL) {
   if (phyloarrayConnection(con)) {
      stm <- paste("SELECT arrayID, barcode, array, sampleID, dye, isDiscarded,",
                   "featureExtraction.* FROM array JOIN hybridisation USING",
                   "(arrayID) JOIN featureExtraction USING (hybridisationID)")
      stm <- paste(stm, expandCondition(condition), sep='')
      rs <- dbSendQuery(con, stm)
      extrs <- fetch(rs, n=-1)
      return(new("extractionList",list(extractions=extrs))) 
   }
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
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



#' choose.samples
#'
#' @param con MySQL connection
#' @param multi multiple selections allowed
#' @param title title
#' @param condition TBA
#' 
#' @return TBA
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

choose.samples <- function (con, multi=TRUE, title='Select samples:', condition=NULL) {
   smps <- fetch.samples(con, condition=condition)
   samples <- select.list(smps$sampleID, multiple=multi, title=title)
   smps <- fetch.samples(con, condition=list(list(field='sampleID',value=samples)))
   return(smps)
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



#' Choosing (and creating) a directory
#' 
#' @param ... parameters to pass
#'
#' @return TBA
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

chooseDir <- function (...) {
  choice = ''
  while (choice == '') {
    choice <- guiDlgDir(dir='',...)
  }
  # if it doesn't exist, it must be made
  if (!file.exists(choice)) {
    dir.create(choice,recursive=TRUE)
  }
  return(choice)
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


#' header part
#' 
#' @param file filename
#' @param title title
#'
#' @return TBA
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

HTMLReportBegin <- function (file="report.html",title="Report Title") {
   cat(paste("<html><head><title>",
   title, "</title></head>",
   "<body>",
   sep = ""), file=file, append=TRUE)
}

#' generic report part
#' 
#' @param x TBA
#' @param file TBA
#' @param CSSfile TBA
#'
#' @return TBA
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

HTMLReport <- function (x, file=paste(file.path(x$directory,x$filename),'html',sep='.'), CSSfile='R2HTML.css') {
   HTMLReportBegin(file)
   file.copy(file.path(system.file(package='R2HTML'),'output','R2HTML.css'),file.path(x$directory,'R2HTML.css'))
   HTMLCSS(file=file, CSSfile=CSSfile)
   HTML(x,file=file)
   HTMLReportEnd(file)
   return(file)
}

#' End html report part
#' 
#' @param file filename
#'
#' @return TBA
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

HTMLReportEnd <- function (file="report.html") {
   cat("<hr size=1></body></html>",
       file=file, append=TRUE)
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities

outlierPvalue <- function (x, avgvar) {
   test <- chisq.out.test(x, variance = avgvar)
   return(test$p.value)
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities


rsd <- function(x,y=NULL) {
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
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
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



