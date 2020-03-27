#' @title Read Simple OTU Tables into a Phyloseq Object
#' @description Read simple OTU tables, mapping and taxonomy files into a
#' \code{\link{phyloseq-class}} object.
#' @details Simple OTU tables, mapping and taxonomy files will be converted
#' to \code{\link{phyloseq-class}}.
#' @param otu.file A simple otu_table with '.csv' extension 
#' @param taxonomy.file A simple taxonomy file with '.csv' extension 
#' @param metadata.file A simple metadata/mapping file with .csv extension
#' @param sep CSV file separator
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples
#' # NOTE: the system.file command reads these example files from the
#' # microbiome R package. To use your own local files, simply write
#' # otu.file <- "/path/to/my/file.csv" etc.
#'
#' #otu.file <-
#' #   system.file("extdata/qiita1629_otu_table.csv",
#' #   package='microbiome')
#'
#' #tax.file <- system.file("extdata/qiita1629_taxonomy_table.csv",
#' #        package='microbiome')
#'
#' #meta.file <- system.file("extdata/qiita1629_mapping_subset.csv",
#' #     package='microbiome')
#'
#' #p0 <- read_csv2phyloseq(
#' #        otu.file=otu.file, 
#' #        taxonomy.file=tax.file, 
#' #        metadata.file=meta.file)
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
read_csv2phyloseq <- function(otu.file=NULL, taxonomy.file=NULL,
    metadata.file=NULL, sep = ",") {

    s.meta <- read.csv(metadata.file, row.names=1, check.names=FALSE, sep = sep)
    s.sampledata <- sample_data(s.meta)
    
    s.otu <- read.csv(otu.file, row.names=1, check.names=FALSE, sep = sep)

    # Ensure that the OTU table is OTU x samples
    if (any(rownames(s.otu) %in% rownames(s.meta))) {
        s.otu <- t(s.otu)
    }

    s.otu.table <- otu_table(s.otu, taxa_are_rows=TRUE)
    s.tax_table <- read_taxtable(taxonomy.file, sep = sep)

    # Add the finest observation level to the tax table
    # (not sufficient to have it as row.names)
    if (!all(rownames(s.tax_table) == s.tax_table[, ncol(s.tax_table)])) {
        s.tax_table <- cbind(s.tax_table, OTU = rownames(s.tax_table))
        s.tax_table <- tax_table(s.tax_table)
    }

    pseq <- merge_phyloseq(s.otu.table, s.tax_table, s.sampledata) 
    return(pseq)
}



