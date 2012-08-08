# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <leo.lahti@iki.fi>. All rights reserved.

# This file is a part of the microbiome R package

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


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


