#' Description: brew embedded documentation helper function
#'
#' Arguments:
#'   @param fname Ouput file name
#'   @param captiontext Caption text
#'   @param width Image width
#'   @param height Image height
#'   @param rotate Optional. Image rotation angle.
#'   @param include.figname Include figure name		     	       
#' Returns:
#'   @return Saves image in the output file 
#'     	      
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

include.fig <- function (fname, captiontext, width = "5cm", height="5cm", rotate = NULL, include.figname = TRUE) {    
  
    cat("\\begin{figure}[h!]\\begin{center}\n")

    if (length(rotate)>0) {cat("\\rotatebox{",rotate,"}{")}
    cat("\\includegraphics")
    cat("[width=",width,", height=",height,"]")
    cat("{")
    cat(fname)
    cat("}")
    if (length(rotate)>0) {cat("}")}
    cat("\n")

    cat("\\caption*{")
    cat(captiontext)
    cat("}")

    cat("\n\\end{center}\\end{figure}\n\n\n")
    
}


#' Description: brew embedded documentation helper function to include 
#' multiple figures
#'
#' Arguments:
#'   @param fnames Ouput file name
#'   @param captiontext Caption text
#'   @param width Image width
#'   @param height Image height
#'   @param rotate Optional. Image rotation angle.
#'		     	       
#' Returns:
#'   @return Saves image in the output file 
#'     	      
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

include.figs <- function (fnames, captiontext, width = "5cm", height="5cm", rotate = NULL) {    

  
    cat("\\begin{figure}[h!]\\begin{center}\n")


    for (i in 1:length(fnames)) {

      fname = fnames[[i]]
      
      if (length(rotate[[i]])>0) {cat("\\rotatebox{",rotate[[i]],"}{")}
      cat("\\includegraphics")
      cat("[width=",width,", height=",height,"]")
      cat("{")
      cat(fname)
      cat("}")
      if (length(rotate)>0) {cat("}")}
      cat("\n")
            
    }
    
    cat("\\caption*{")
    cat(captiontext)
    cat("}")

    cat("\n\\end{center}\\end{figure}\n\n\n")
    
}
