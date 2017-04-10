#' @title  Saving Various Files as Comma Seperated Files
#' @description This command will save the otu,taxonomy and meta data as *.csv file for manual inspection.
#' @param x \code{\link{phyloseq-class}} object
#' @param type "OTU" or "TAXA" or "METADATA"
#' @return  Output file path (a string)
#' @export
#' @examples \dontrun{
#'   library(microbiome)
#'   data(dietswap)
#'   pseq <- dietswap
#'   write_phyloseq(pseq, "OTU")
#'   write_phyloseq(pseq, "TAXA")
#'   write_phyloseq(pseq, "METADATA")
#'  }
#' @keywords utilities
write_phyloseq <- function(x, type){
  if (type == "OTU"){
    f <- paste(getwd(), "otu_table.csv", sep = "/");
    message("Writing OTU in the file ", f);
    y <- as.data.frame(x@otu_table);
    write.csv(y, file = f, fileEncoding = "UTF-16LE");
    } else if (type == "TAXA"){
      f <- paste(getwd(), "taxa_table.csv", sep = "/");
      message("Writing TAXA in the file ", f)
    y <- as.data.frame(x@tax_table);
    write.csv(y, file = f, fileEncoding = "UTF-16LE");
    } else if (type == "METADATA"){
      f <- paste(getwd(), "meta_table.csv", sep = "/");
      message("Writing METADATA in the file ", f)
    y <- as.data.frame(x@sam_data);
    write.csv(y, file = f, fileEncoding = "UTF-16LE")
    }
  
}
