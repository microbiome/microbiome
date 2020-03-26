#' @title Default Colors
#' @description Default colors for different variables.
#' @param x Name of the variable type ("Phylum")
#' @param v Optional. Vector of elements to color.
#' @return Named character vector of default colors
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("microbiome")
#' @export
#' @examples col <- default_colors("Phylum")
#' @keywords utilities
default_colors <- function (x, v=NULL) {

    if (x == "Phylum") {
        #http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
        #https://www.r-graph-gallery.com/42-colors-names/
        col <- c(
            "Actinobacteria" = "darkgreen",
            "Firmicutes" = "blue",
            "Proteobacteria" = "black",
            "Verrucomicrobia" = "darkblue",
            "Bacteroidetes" = "red",
            "Spirochaetes" = "darkgray",
            "Fusobacteria" = "lightblue",
            "Cyanobacteria" = "deepskyblue3")
    }

    if (x == "Genus") {
        col <- c(
        "Faecalibacterium" = "blue",
        "Bifidobacterium" = "magenta",
        "Escherichia" = "black",
        "Proteobacteria" = "black",
        "Streptococcus" = "black",
        "Enterococcus" = "black",
        "Staphylococcus" = "black",
        "Pathogens" = "black",    
        "Alistipes" = "lightblue",    
        "Other" = "darkgray",
        "Prevotella" = "darkgreen",
        "Bacteroides" = "red",
        "Ruminococcus" = "darkblue",
        "Oscillospira" = "pink"
        )
    }  

    col

}

