#' @title estimate_variability
#' @description Variability analysis. 
#'
#' @details 
#' Average correlation between samples in the input data within each
#' group with the overall group-wise average. Picks the lower triangular
#' matrix to avoid duplicating the correlations. Returns correlations and
#' stability estimate (inter-individual; average of the
#' correlations). Can also be used to calculate temporal variability
#' between two data sets (intra-individual), given appropriate sample
#' metadata.
#'
#' @param x phyloseq object with the OTU matrix and sample metadata. The sample metadata should contain the fields "sample" and "group" (or another grouping variable specified in the group_by argument) for the "interindividual" method. For the "intraindividual" method, also provide "time", and "subject" fields.
#' @param type Variability type: 'interindividual' or 'intraindividual'
#' @param group_by variable to be used in grouping. By default: "group"
#' @param method correlation method (see ?cor)
#'
#' @return List with correlations, group-wise statistics, and ANOVA 
#'         linear model p-value for group differences.
#'
#' @export
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
#' @examples 
#' \dontrun{
#' # Example data
#' library(microbiome)
#' x <- download_microbiome("dietswap")
#' sample_data(x)$time <- sample_data(x)$timepoint.within.group
#' # Estimate inter-individual variability
#' res <- estimate_variability(x, "interindividual")
#' # Estimate intra-individual variability
#' res <- estimate_variability(x, "intraindividual")
#' }
#' @references 
#' The inter- and intra-individual variability are calculated
#' as described in Salonen et al. ISME J. 8:2218-30, 2014.
#' To cite this R package, see citation('microbiome')
#' 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
estimate_variability <- function(x, type, group_by = "group", method = "spearman") {

    # Pick metadata		     
    meta <- sample_data(x)
    
    # OTU data Log10
    otu <- log10(t(otu_table(x)@.Data))

    correlation <- NULL
 		         
    # Split the data by group
    group <- NULL
    if (!group_by %in% names(meta)) {
      meta[[group_by]] <- rep("completedata", nrow(meta))
    }
    datasets <- split(as.data.frame(otu), meta[[group_by]], drop = TRUE)

    if (type == "interindividual") {
    
      tmp <- setdiff(c("sample", group_by), names(meta))
      if (length(tmp) > 0) {
        stop(paste("The following variables needed by estimate_variability function type=interindividual are missing from sample metadata:", paste(tmp, collapse = ",")))
      }
  
      # Within-matrix stability NOTE: earlier this was calculated as
      # the average of upper triangular correlation matrix This is
      # heavily biased since the values are dependent Now replaced
      # by calculating correlations against the mean of the whole
      # sample set cors <- lower.triangle(cor(dat1))
      dfs <- NULL      
      for (ds in names(datasets)) {
        dat1 <- datasets[[ds]]
        cors <- as.vector(cor(t(dat1), matrix(colMeans(dat1)), method = method, use = "pairwise.complete.obs"))
        dfs <- rbind(dfs, data.frame(group = rep(ds, length(cors)),
                   	             sample = rownames(dat1),
  		   		     correlation = cors))
	#variability[[ds]] <- list(correlations = cors, variability = mean(cors))
      }
  
      pval <- anova(lm(correlation ~ group, data = dfs))[["Pr(>F)"]][[1]]
      stats <- dfs %>% group_by(group) %>% summarize(mean = mean(correlation), sd = sd(correlation))
      variability <- list(data = dfs, statistics = stats, p.value = pval)

    } else if (type == "intraindividual") {

      tmp <- setdiff(c("time", "subject", "sample", group_by), names(meta))
      if (length(tmp) > 0) {
        stop(paste("The following variables needed by estimate_variability function type=intraindividual are missing from sample metadata:", paste(tmp, collapse = ",")))
      }

      if (!all(sapply(split(meta$time, meta[[group_by]]), function (x) {length(unique(x))}) == 2)) {
        stop("Two time points needed for each group for the intraindividual type. Some groups are having a different number of time points.")
      }

      variability <- list()
      dfs <- NULL
      for (ds in names(datasets)) {
      
        # Pick the data and metadata for this group
        xsub <- datasets[[ds]]        
	msub <- meta[rownames(xsub),]

        # Use interindividual functionality to assess correlations
	# within subjects. Subjects are used as groups
        datasets2 <- split(xsub, droplevels(msub[["subject"]]))
	cors <- c()
        for (subj in names(datasets2)) {
          dats <- datasets2[[subj]]
          cors[[subj]] <- cor(unlist(dats[1,]), unlist(dats[2,]), method = method, use = "pairwise.complete.obs")
        }

        dfs <- rbind(dfs, data.frame(group = rep(ds, length(cors)),
                   	             subject = names(cors),
  		   		     correlation = cors))

      }

      # Between time point correlations within subjects
      # and the mean over those correlations
      pval <- anova(lm(correlation ~ group, data = dfs))[["Pr(>F)"]][[1]]
      stats <- dfs %>% group_by(group) %>% summarize(mean = mean(correlation, na.rm = T), sd = sd(correlation, na.rm = T))
      variability <- list(data = dfs, statistics = stats, p.value = pval)

    }
    
    variability
    
}


