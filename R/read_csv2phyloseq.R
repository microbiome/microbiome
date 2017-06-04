#' @title Read Simple OTU Tables into a Phyloseq Object
#' @description Read simple OTU tables, mapping and taxonomy files into a
#' \code{\link{phyloseq-class}} object.
#' @details Simple OTU tables, mapping and taxonomy files will be converted
#' to \code{\link{phyloseq-class}}.
#' @param otu.file A simple otu_table with '.csv' extension 
#' @param taxonomy.file A simple taxonomy file with '.csv' extension 
#' @param metadata.file A simple metadata/mapping file with .csv extension
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples
#' \dontrun{
#' # TODO: add example files in inst/extdata/
#' # and read here with
#' # system.file('inst/extdata/..', package='microbiome') command
#' # To make this example executable
#' p0 <- read_csv2phyloseq(otu.file='otu_table.csv', 
#'     taxonomy.file='taxa_table.csv', 
#'     metadata.file='meta_table.csv')
#' }
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
read_csv2phyloseq <- function(otu.file=NULL, taxonomy.file=NULL,
    metadata.file=NULL) {
    s.otu <- read.csv(otu.file, row.names=1, check.names=FALSE)
    s.otu.table <- otu_table(s.otu, taxa_are_rows=TRUE)
    s.tax <- read.csv(taxonomy.file, row.names=1, check.names=FALSE)
    s.taxmat <- as.matrix(s.tax)
    s.tax_table=tax_table(s.taxmat)
    s.meta <- read.csv(metadata.file, row.names=1, check.names=FALSE)
    s.smapledata=sample_data(s.meta)
    simple_pseq <- merge_phyloseq(s.otu.table, s.tax_table, s.smapledata)    
    return(simple_pseq)
}

