# Copyright (C) 2006-2012 Douwe Molenaar, Janne Nikkila, 
# Leo Lahti, and Jarkko Salojarvi 
#
# Contact: <leo.lahti@iki.fi>. All rights reserved.
#
# This file is a part of the microbiome R package
#
# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' Impute missing values from a Gaussian. 
#' 
#' @param X data matrix (features x samples)
#'
#' @return imputed data matrix
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples #X2 <- impute(X) 
#' @keywords utilities

impute <- function (X) {

  for (i in 1:ncol(X)) {
    x <- X[, i]
    nas <- is.na(x)
    X[nas, i] <- rnorm(sum(nas), mean(x[!is.na(x)]), sd(x[!is.na(x)]))
  }
  
  X
  
}


#' Strip string i.e. remove spaces from the beginning and end
#' @param s string or character vector
#'
#' @return Stripped string
#' @export 
#' @references
#' See citation("microbiome")
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#' #s2 <- Strip(s) 
#' @keywords utilities
Strip <- function (s) {

  ss <- c()
  
  for (i in 1:length(s)) {

    si <- s[[i]]
    if (!is.na(si)) {
      # Strip string i.e. remove spaces from the beginning and end
      while (substr(si,1,1)==" ") {
        si <- substr(si, 2, nchar(si))
      }
      while (substr(si, nchar(si), nchar(si))==" ") {
        si <- substr(si, 1, nchar(si) - 1)
      }
    }
    ss[[i]] <- si
  }

  ss
}




#' Description: Sort data frame dd by columns like: esort(dd, -z, b)
#'
#' Arguments:
#'   @param x data frame to sort
#'   @param sortvar sorted variable/s
#'   @param ... further parameters to pass
#'
#' Returns:
#'   @return sorted data frame
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

esort <- function(x, sortvar, ...) {

  attach(x)
  x <- x[with(x,order(sortvar,...)),]
  return(x)
  detach(x)
}


#' Description: This feature for dendrograms took like an hour to find out
#'               because of unadequate documentation...
#'
#' Arguments:
#'   @param dd TBA
#'
#' Returns:
#'   @return dd
#'
#' @export 
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.dd.order <- function(dd){
  unlist(dd)
}

#' Description: R or S function to write an R list to an ASCII file.
#'   This can be used to create files for those who want to use
#'   a spreadsheet or other program on the data.
#'
#' Arguments:
#'   @param x TBA
#'   @param file TBA
#'   @param append TBA
#'
#' Returns:
#'   @return TBA
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @note Modified from version MHP July 7, 2004. 
#' @keywords utilities

list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep=""), append = F) {

   tmp.wid = getOption("width")  # save current width
   options(width=10000)          # increase output width
   sink(file, append=append)     # redirect output to file
   print(x)                      # print the object
   sink()                        # cancel redirection
   options(width=tmp.wid)        # restore linewidth
   return(invisible(NULL))       # return (nothing) from function
}



#' Description: 
#' Get lower triangle of a square matrix 
#' as a numeric vector such that
#' row-by-row, picking elements in the order
#' 2,1;3,1;3,2;4,1,...
#'        
#' Arguments:
#'   @param mat data matrix
#' Returns:
#'   @return lower triangle as vector 
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
lower.triangle <- function (mat) {

        elements<-c()
        nr=dim(mat)[[1]]
        nc=dim(mat)[[2]]
        for (i in 2:nr) {
                for (j in 1:(i-1)) {
                        elements<-c(elements,mat[i,j])
                }
        }       
        elements
}
