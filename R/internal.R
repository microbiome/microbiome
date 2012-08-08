# database utilities for package-internal use only

#' Tests whether the database connection is a phyloarray connection.
#' Expands one element (one "field", "value" pair list) from a list 
#' of "field", "value" pair lists
#'
#' @param con a MySQL database connection.
#'
#' @return TRUE when the test succeeds. Otherwise a program halt.
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

phyloarrayConnection <- function (con) {
   if (!(class(con)=='MySQLConnection')) {
      stop('Input must be a DBI connection to a phyloarray database')
   }
   essential <- c("array",
                  "arraydesign",
                  "arrayfeature",
                  "arrayhybridisations",
                  "featureextraction",
                  "featuremeasurement",
                  "hybridisation",
                  "image",
                  "oligo",
                  "oligoclass",
                  "oligotargetpair",
                  "phylogeny",
                  "project",
                  "sample",
                  "slide",
                  "target",
                  "taxon",
                  "taxonlevel",
                  "taxtotax")
   if (!(length(intersect(dbListTables(con),essential))==length(essential))) {
      stop('Essential tables missing in the connected database. Not a phyloarray database?')
   }
   return(TRUE)
}


#' Expands one element (one "field", "value" pair list) from a list of "field", "value" pair lists
#' @param elm TBA
#'
#' @return TBA
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

expandElement <- function (elm) {
   if (is.list(elm)) {
      if (is.null(elm$field)) {
         stop("Missing 'field' field in condition")
      }
      if (!is.vector(elm$field)) {
         stop("Field 'field' must be vector with length 1")
      }
      if (length(elm$field)>1) {
         stop("Field 'field' must be vector with length 1")
      }
      if (!is.vector(elm$value)) {
         stop("'value' field must be vector")
      }
      if (is.null(elm$value)) {
         stop("Missing 'value' field in condition")
      }
      if (is.character(elm$value)) {
         elm$value <- paste("'",elm$value,"'",sep='')
      }
       return(paste("(",paste(elm$field,elm$value,sep='=',collapse=' OR '),")",sep=''))
   }
   else {
      stop("Argument 'condition' must be a list of lists with 'field' and 'value' pairs")
   }
}  

#' Expands a list of "field", "value" pair lists into an SQL condition
#' 
#' Arguments:
#'  @param condition condition
#'
#' Returns:
#'  @return TBA
#'
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # TBA
#' @keywords utilities


expandCondition <- function (condition) {
   if (is.null(condition)) {
      return(NULL)
   } else {
      if (is.list(condition)) {
         return(paste(" WHERE",paste(lapply(condition, expandElement),collapse=' AND ')))
      } else {
         stop("Argument 'condition' must be a list of lists with 'field' and 'value' pairs")
      }
   }  
}


#' Description: populate radiobuttons 
#'
#' Arguments:
#'   @param tt TBA
#'   @param title TBA
#'   @param var.names TBA 
#'   @param var.values TBA
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

populate.radiobuttons <- function(tt, title, var.names, var.values, var.init) {

  title.font <- tkfont.create(weight = "bold", size = 10)

  frm <- tkframe(tt, relief = "groove", borderwidth = 2)

  label.widget <- tklabel(frm, text = title,font = title.font)

  tkpack(label.widget, side = "top")

  for (i in 1:length(var.values)){
    button.widget <- tkradiobutton(frm, text = var.names[i], 
                                  variable = var.init, value = var.values[i])
    tkpack(button.widget,side = "top")
  }

  tkpack(frm,side = "top")

  return(list(frame = frm, var = var.init))

}


