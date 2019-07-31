#' @title Download microbiome data sets
#' @description Download microbiome data sets
#' @param id Data set name. For options, see download_microbiome()
#' @param format Data format ("phyloseq" or "original")
#' @return Data set
#' @examples # x <- download_microbiome("peerj32")
#' @export
#' @references 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
download_microbiome <- function (id = "datasets", format = NULL) {

  if (id == "datasets") {
    datasets <- c("atlas1006", "dietswap", "peerj32")
    message(paste("The available data sets ids are: ", paste(datasets, collapse = ", ")))
    return(datasets)
  }

  if (id == "atlas1006") {

    data <- download_atlas(format)

  } else if (id == "dietswap") {

    data <- download_dietswap()
  
  } else if (id == "peerj32") {

    data <- download_peerj32()

  }

  data

}


download_peerj32 <- function (...) {

    message("Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/")

    peerj32 <- list()
    data.dir <- system.file("inst/extdata", package = "microbiome")
    for (nam in c("lipids", "microbes", "meta")) {
      peerj32[[nam]] <- read.table(paste0(data.dir, "/peerj32_", nam, ".csv"), 
         sep = "\t", header = TRUE, row.names = 1)
    }

    data <- peerj32

    # The formatting of taxon names needs to be harmonized for microbiome package:
    colnames(data$microbes) <- gsub("\\.", " ", colnames(data$microbes))
    # Put back 'rel.' periods
    colnames(data$microbes) <- gsub("rel $", "rel.", colnames(data$microbes))
    colnames(data$microbes) <- gsub("Clostridium  sensu stricto ", "Clostridium (sensu stricto)", colnames(data$microbes))
    # Convert to matrix 
    data$microbes <- as.matrix(data$microbes)

    # Harmonize the field names etc.
    data$meta <- harmonize_fieldnames(data$meta)

    # Harmonize field contents
    data$meta <- suppressWarnings(harmonize_fields(data$meta))

    # Convert in phyloseq format
    physeq <- hitchip2physeq(data$microbes, data$meta)

    data$phyloseq <- physeq
    
    return(data)     #list(physeq = physeq, data = data)

}





#' @title Download HITChip Atlas
#' @description Download HITChip Atlas 
#' @param ... Arguments to pass
#' @return Data set
#' @examples \dontrun{download_atlas()}
#' @importFrom rdryad download_url
#' @references 
#'   Lahti et al. Tipping elements of the human intestinal ecosystem. 
#'   Nature Communications 5:4344, 2014.
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
download_atlas <- function (...) {

  message("Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from 
  Data Dryad: http://doi.org/10.5061/dryad.pk75d")

  # Define the data URL
  url <- download_url('10255/dryad.64665')[[1]]

  # Download the data
  data <- read.table(url, sep = "\t", row.names = 1, header = TRUE)

  # The data comes from public repository (Dryad). 
  # The formatting of taxon names needs to be harmonized for microbiome package:
  # Fix some broken names from the original release..
  # ie. replace 'Clostridium..sensu.stricto.les' with 'Clostridiales'
  colnames(data) <- gsub("Clostridium..sensu.stricto.les", "Clostridiales", colnames(data))
  # Remove periods
  colnames(data) <- gsub("\\.", " ", colnames(data))
  # Put back 'rel.' periods
  colnames(data) <- gsub("rel $", "rel.", colnames(data))
  colnames(data) <- gsub("Clostridium  sensu stricto ", "Clostridium (sensu stricto)", colnames(data))

  # Convert to matrix 
  otu <- as.matrix(data)

  url <- download_url('10255/dryad.64666')[[1]]
  meta <- unique(read.table(url, sep = "\t", row.names = 1, header = TRUE))

  # Add SampleIDs as a separate column from rownames
  meta$sample <- rownames(meta)

  # Order BMI groups in correct order
  # (see README at http://datadryad.org/resource/doi:10.5061/dryad.pk75d 
  # for details)

  # Harmonize the field names etc.
  colnames(meta) <- harmonize_fieldnames(colnames(meta))

  # Harmonize field contents
  meta <- harmonize_fields(meta)

  # Convert in phyloseq format
  physeq <- hitchip2physeq(otu, meta)
  res <- physeq

  # Remove a single sample per subject/timepoint combination
  # since this R package version is mainly for pedagogical purposes
  # and the original data is available at the original publication
  res <- collapse_replicates(res, replicate_fields = c("subject", "time"))
  
  res

}



#' Download Diet Swap Study
#'
#' @param ... Arguments to be passed
#' @return Data set
#' @examples \dontrun{download_dietswap()}
#' @importFrom rdryad download_url
#' @references 
#'   O'Keefe et al. Nature Communications 6:6342, 2015.
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
download_dietswap <- function (...) {

  message("Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n")

  # Define the data URL
  url <- unique(download_url('10255/dryad.78878'))

  # Download the data
  data <- read.table(url, sep = ",", row.names = 1, header = TRUE)

  # The data comes from public repository (Dryad). 
  # The formatting of taxon names needs to be harmonized for microbiome package:
  # Fix some broken names from the original release..
  # ie. replace 'Clostridium..sensu.stricto.les' with 'Clostridiales'
  colnames(data) <- gsub("^Clostridia$", "Clostridium (sensu stricto)", colnames(data))
  # Remove periods
  colnames(data) <- gsub("\\.", " ", colnames(data))
  # Put back 'rel.' periods
  colnames(data) <- gsub("rel $", "rel.", colnames(data))

  # Convert to matrix 
  data <- as.matrix(data)

  # -------------------------------------------

  url <- unique(download_url('10255/dryad.78880'))
  meta <- read.table(url, sep = ",", row.names = 1, header = TRUE)

  # Add SampleIDs as a separate column from rownames
  meta$sample <- rownames(meta)

  # Order BMI groups in correct order
  # (see README at http://datadryad.org/resource/doi:10.5061/dryad.pk75d 
  # for details)

  # Harmonize time
  meta$timepoint <- meta$timepoint.total
  meta$timepoint.total <- NULL

  meta$timepoint.within.group <- meta$timepoint.group
  meta$timepoint.group <- NULL

  # Harmonize the field names etc.
  # colnames(meta) <- suppressWarnings(harmonize_fieldnames(colnames(meta)))

  # Harmonize field contents
  meta <- suppressWarnings(harmonize_fields(meta))

  # Collect the atlas data and metadata into a single object
  # atlas <- list(microbes = data, meta = meta)

  # Convert in phyloseq format
  pseq <- hitchip2physeq(data, meta)

  pseq

}
