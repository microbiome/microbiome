#' @title ANOVA for phyloseq
#' @description ANOVA test for multi-group comparison for phyloseq objects.
#' @param x \code{\link{phyloseq-class}} object 
#' @param group Metadata field specifying the groups.
#' @param p.adjust.method p-value correction method for p.adjust function 
#'               (default 'BH'). For other options, see ?p.adjust
#' @return Corrected ANOVA p-values for multi-group comparison. Uncorrected post-hoc p-values for each pairwise comparison. Group-wise averages.
#' @examples 
#'   data(peerj32)
#'   pval <- check_anova(peerj32$phyloseq, "group")
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_anova <- function (x, group, p.adjust.method = "BH") {

  . <- p.anova <- NULL

  # We need taxa x samples matrix
  mydata <- get_taxa(x)
  if (!taxa_are_rows(x)) {mydata <- t(mydata)}

  metadata <- sample_data(x)
  metadata$group <- metadata[[group]]

  # For each taxon, check if group
  
  # ANOVA p-values
  pv <- matrix(NA, nrow = nrow(mydata), ncol = 1 + choose(length(levels(metadata$group)), 2))
  rownames(pv) <- rownames(mydata)

  for (tax in rownames(mydata)) {

    dfs <- metadata
    dfs$signal <- mydata[tax, rownames(dfs)]

    # fit  <- anova(lm(signal ~ group, data = dfs))
    # Type I sequential SS
    fit <- aov(terms(signal ~ group), data = dfs)     
    pv1 <- summary(fit)[[1]]["group", "Pr(>F)"]

    # Group differences, confidence intervals and p-values
    tukey.table <- TukeyHSD(fit)$group
    pv2 <- tukey.table[, "p adj"]
    names(pv2) <- paste("p.", names(pv2), sep = "")

    # Type III results compare each term with the full model.
    # Alternatively, we can use anova(fit.model1, fit.model2) to compare nested models directly.     
    # pv[[tax]] <- drop1(fit,~.,test="F")["group", "Pr(>F)"]

    # ANOVA all groups p-value and then the pairwise posthocs
    pv[tax,] <- c(pv1, pv2)
    
  }
  colnames(pv) <- c("p.anova", names(pv2))

  # Adjust ANOVA group-level p-values
  padj <- p.adjust(pv[, "p.anova"], method = p.adjust.method)  
  taxa <- names(padj)

  # Signal group-level means
  dfm <- as.data.frame(t(mydata))
  dfm$group <- metadata$group

  # With dplyr
  relmeans <- dfm %>% group_by(group) %>% summarize_each(funs(mean(na.omit(.))))
  rnams <- as.character(relmeans$group)
  relmeans <- as.matrix(relmeans)[, -match("group", colnames(relmeans))]
  relmeans <- t(apply(relmeans, 2, as.numeric))
  colnames(relmeans) <- paste("ave.", rnams, sep = "")    

  # Combine the data frames
  res <- as.data.frame(cbind(taxa = taxa,
                             padj = padj,
			     relmeans))

  res$padj <- as.numeric(as.character(res$padj))
  for (i in 3:ncol(res)) {
    res[,i] <- as.numeric(as.character(res[,i]))
  }

  
  # Order by p 
  # res <- esort(res, p.anova)
  res <- res %>% arrange(padj)
  
  res

}


