#' Description: Default list of removed phylotypes and oligos
#'
#' Arguments:
#'  @param chip Chip name (HIT/MIT/PIT/Chick)Chip
#' Returns:
#'   @return List of removed oligos and phylotypes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

phylotype.rm.list <- function (chip) {

  rm.phylotypes <- list()

  if (chip == "HITChip") {
    
    rm.phylotypes[["oligos"]] <- c("UNI 515", "HIT 5658", "HIT 1503", "HIT 1505", "HIT 1506")
    rm.phylotypes[["species"]] <- c("Victivallis vadensis")
    rm.phylotypes[["level 1"]] <- c("Lentisphaerae")
    rm.phylotypes[["level 2"]] <- c("Victivallis")

  } else if (chip == "MITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["level 1"]] <- c()
    rm.phylotypes[["level 2"]] <- c()

  } else if (chip == "PITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["level 1"]] <- c()
    rm.phylotypes[["level 2"]] <- c()

  } else if (chip == "ChickChip") {
    warning("No universal probes excluded from ChichChip yet!")
  }

  rm.phylotypes

}





#' Description: List scaling methods
#'
#' Arguments:
#'
#' Returns:
#'   @return List of scaling methods
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.scaling.methods <- function () {

  list('none'='none',
                #'minimum/median'='minmed',
                'minimum/maximum'='minmax',
                'minmax'='minmax',
                #'median'='med',
                'quantile'='quant'
                #'normExp+MedianFC'='normExpMedianFC',
                #'normExp+quantile'='normExpQuant'
   )

}


#' Description: List clustering metrics
#'
#' Arguments:
#'
#' Returns:
#'   @return list of clustering metrics
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.clustering.metrics <- function () {

  list('Pearsons correlation coefficient'='correlation',
                 'euclidian'='euclidian')
}



#' Description: List color scales
#'
#' Arguments:
#'
#' Returns:
#'   @return list of color scales
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.color.scales <- function () {
  ## Different colour scales
  list('white/blue'=colorRampPalette(c("white","darkblue"),interpolate='linear')(100),
       'white/black'=colorRampPalette(c("white","black"),interpolate='linear')(100),
       'black/yellow/white'=colorRampPalette(c("black","yellow","white"),bias=0.5,interpolate='linear')(100))

}


