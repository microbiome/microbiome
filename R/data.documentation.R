#' @title Probiotics intervention data documentation 
#' @description The peerj32 data set contains high-through profiling data from 
#' 389 human blood serum lipids and 130 intestinal genus-level bacteria from 
#' 44 samples (22 subjects from 2 time points; before and after 
#' probiotic/placebo intervention). The data set can be used to investigate 
#' associations between intestinal bacteria and host lipid metabolism. 
#' For details, see \url{http://dx.doi.org/10.7717/peerj.32}.
#' @name peerj32
#' @docType data
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com} 
#' @references Lahti et al. (2013) PeerJ 1:e32 \url{http://dx.doi.org/10.7717/peerj.32}
#' @usage data(peerj32)
#' @format List of the following data matrices as described in detail in Lahti et al. (2013):
#' \itemize{
#'   \item lipids: Quantification of 389 blood serum lipids across 44 samples
#'   \item microbes: Quantification of 130 genus-like taxa across 44 samples
#'   \item meta: Sample metadata including time point, gender, subjectID, sampleID and treatment group (probiotic LGG / Placebo)
#'   \item phyloseq The microbiome data set converted into a \code{\link{phyloseq-class}} object.
#' }
#'   
#' @keywords data
NULL 


#' @title Diet swap study data 
#' @description The diet swap data set represents a study with African and African American groups undergoing a two-week diet swap.
#' For details, see \url{http://www.nature.com/ncomms/2015/150428/ncomms7342/full/ncomms7342.html}.
#' @name dietswap
#' @details The data is also available for download from the Data Dryad repository \url{http://datadryad.org/resource/doi:10.5061/dryad.1mn1n}. 
#' @docType data
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com} 
#' @references 
#'   O'Keefe et al. Nature Communications 6:6342, 2015.
#'   \url{http://www.nature.com/ncomms/2015/150428/ncomms7342/full/ncomms7342.html}
#'   To cite the microbiome R package, see citation('microbiome') 
#' @usage data(dietswap) 
#' @format The data set in \code{\link{phyloseq-class}} format.  
#' @keywords data
NULL 



#' @title HITChip Atlas data (n = 1006)
#' @description This data set contains genus-level microbiota profiling with HITChip for 1006 western adults with no reported health complications, reported in Lahti et al. (2014) \url{http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html}.
#' @name atlas1006
#' @details The data is also available for download from the Data Dryad repository \url{http://doi.org/10.5061/dryad.pk75d}. 
#' @docType data
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com} 
#' @references 
#'   Lahti et al. Tipping elements of the human intestinal ecosystem. 
#'   Nature Communications 5:4344, 2014.
#'   To cite the microbiome R package, see citation('microbiome') 
#' @usage data(atlas1006) 
#' @format The data set in \code{\link{phyloseq-class}} format.  
#' @keywords data
NULL 


#' @title HITChip taxonomy
#' @description HITChip taxonomy table
#' @name hitchip.taxonomy
#' @docType data
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com} 
#' @references 
#'   Lahti et al. Tipping elements of the human intestinal ecosystem. 
#'   Nature Communications 5:4344, 2014.
#'   To cite the microbiome R package, see citation('microbiome') 
#' @usage data(hitchip.taxonomy)
#' @format List with the element 'filtered', including a simplified version of the HITChip taxonomy.
#' @keywords data
NULL 