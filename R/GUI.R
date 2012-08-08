# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <leo.lahti@iki.fi>. All rights reserved.

# This file is a part of the microbiome R package

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: low.level functions required by GUI
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

