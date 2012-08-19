#' Description: Computes the given lme model for each variable (row) in the given data frame
#'      Designed for HITchip matrices, but is applicable to any other matrix also.
#'      NOTE: does not take log automatically!
#'
#' Arguments:
#'   @param fixed TBA
#'   @param vars TBA
#'   @param d.matrix TBA
#'   @param d.info TBA
#'   @param random TBA
#'   @param correlation TBA
#'   @param weights TBA
#'   @param subset TBA
#'   @param method REML or ML
#'   @param na.action TBA
#'   @param control TBA
#'   @param contrasts TBA
#'   @param keep.data TBA
#'
#' Returns:
#'   @return lmeMatrixRes
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

lmeMatrix <- function(fixed, vars, d.matrix, d.info, random,
                        correlation = NULL, weights = NULL, subset, 
			method = c("REML", "ML"),                                                                      
                        na.action = na.fail, control = list(), contrasts = NULL,
                        keep.data = TRUE) {

  ## initialize
  vars <- rownames(d.matrix)
  lmeMatrixRes <- list()

  message("Computing lm for the following variables and data:\n")
  print(vars)
  print(str(d.info))
  print(str(t(d.matrix)))
  
  for(v in vars){
    print(v)
    data <- d.info
    data$y <- t(d.matrix[v,])

    lmeMatrixRes$lme[[v]] <- try(lme(fixed, data, random, correlation, weights, subset,
                         method, na.action, control, contrasts, keep.data))
    lmeMatrixRes$anova[[v]] <- try(anova(lmeMatrixRes$lme[[v]]))
    lmeMatrixRes$summary[[v]] <- try(summary(lmeMatrixRes$lme[[v]]))
    tmp <- try(intervals(lmeMatrixRes$lme[[v]]))

    lmeMatrixRes$intervals[[v]] <- as.list(tmp)
    
  }
  lmeMatrixRes
}



#' Description:  Computes the given lm model for each variable (row) in the 
#'               given matrix/data frame. Designed for HITchip matrices, 
#' 		 but is applicable to any other matrix also.
#'  		 NOTE: does not take log! 
#'
#' FIXME: merge with lm.matrix2
#'
#' Arguments:
#'   @param formula formula
#'   @param d.matrix data matrix
#'   @param d.info information 
#'   @param type information to return
#'
#' Returns:
#'   @return matrix
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

lm.matrix <- function(formula, d.matrix, d.info, type = "coef") {

  require(limma)

  lm.res <- NULL

  for (v in colnames(d.matrix)) {

      df <- data.frame(list(d.info, y = as.vector(d.matrix[,v])))
      lmfit <- lm(formula, df)
      if (type == "coef") {
        lm.res <- rbind(lm.res, lmfit$coefficients)
      } else if (type == "qval") {
        lm.res <- c(lm.res, anova(lmfit)[["Pr(>F)"]][[1]])
      }
  }

  if (type == "coef") {
    rownames(lm.res) <- colnames(d.matrix)
  } else if (type == "qval") {
    names(lm.res) <- colnames(d.matrix)
    require(qvalue)
    lm.res[is.na(lm.res)] <- 1 # pvalue = 1 for missing vals
    lm.res <- qvalue(lm.res, pi0.method = "bootstrap")$qvalue
  }

  lm.res

}



#' Description: Get Pvals From Lm Matrix
#'
#' Arguments:
#'   @param lmMatrix TBA
#'
#' Returns:
#'   @return anova.pvals
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

getPvalsFromLmMatrix <- function(lmMatrix){
  getP <- function(x){
    if(!is.element("p-value",names(x))) {
      tmp <- t(x['Pr(>F)'])
      tmp <- tmp[1:(length(tmp)-1)]
      names(tmp) <- rownames(x)[1:length(tmp)]
      return(tmp)
    } else {
      if(is.element("p-value",names(x))) {
        tmp <- t(x['p-value'])
        tmp <- tmp[2:length(tmp)]
        names(tmp) <- rownames(x)[2:(length(tmp)+1)]
        return(tmp)
      } else {
        print("Class of supposed lm/lme-object not known!")
      }
    }
  }
  anova.pvals <- sapply(lmMatrix[["anova"]], getP)
  ##rownames(anova.pvals) <- rownames(lmMatrix[[2]][[1]]['Pr(>F)'])
  return(anova.pvals)
}


