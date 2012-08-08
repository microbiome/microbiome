# Copyright (C) 2006-2012 Douwe Molenaar, Janne Nikkilä, Leo Lahti, and 
# Jarkko Salojarvi 
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


#' Fetch projects from the phyloarray database
#' @param con MySQL connection
#' @param condition list of lists with field-value pairs
#'
#' Fetch complete records from the \emph{projects}
#' table. If \code{condition=NULL}; will fetch all records and if
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
  choice <- ''
  while (choice == '') {
    choice <- guiDlgDir(dir = '', ...)
  }
  # if it doesn't exist, it must be made
  if (!file.exists(choice)) {
    dir.create(choice, recursive = TRUE)
  }
  return(choice)
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





