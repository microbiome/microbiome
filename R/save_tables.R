#' @title  Saving Various Files as Comma Seperated Files
#' @description This command will save the otu,taxonomy and meta data as *.csv file for manual inspection.
#' @param x \code{\link{phyloseq-class}} object
#' @param type "OTU" or "TAXA" or METADATA
#' @return  Output file path (a string)
#' @export
#' @examples \dontrun{
#'   data(DynamicsIBD)
#'   p0 <- DynamicsIBD
#'   p0.f <- format_phyloseq(p0)
#'   save_tables(p0.f, "OTU")
#'   save_tables(p0.f, "TAXA")
#'   save_tables(p0.f, "METADATA")
#'  }
#' @keywords utilities
save_tables <- function(x, type){
  if (type == "OTU"){
    message("Writing file OTU file in ", getwd())
    y <- as.data.frame(x@otu_table);
    write.csv(y, file = "otu_table.csv", fileEncoding = "UTF-16LE");
    } else if (type == "TAXA"){
      message("Writing file TAXA file in ", getwd())
    y <- as.data.frame(x@tax_table);
    write.csv(y, file ="taxa_table.csv", fileEncoding = "UTF-16LE");
    } else if (type == "METADATA"){
      message("Writing file METADATA file in ", getwd())
    y <- as.data.frame(x@sam_data);
    write.csv(y, file = "meta_table.csv", fileEncoding = "UTF-16LE")
    }
  paste0("File saved in folder ", getwd())
}
