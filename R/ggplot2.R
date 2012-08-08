#' Description: Barplot of top findings from pairwise.comparisons function
#'              
#' Arguments:
#'   @param top.findings Ouput from pairwise.comparisons function
#'   @param topN number of top findings to plot
#'   @param annot annotation matrix    
#'   @param smat species abundancy matrix
#' Returns:
#'   @return List with top findings from pairwise comparisons and their q-values
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

top.barplots <- function (top.findings, topN = 5, annot, smat) {

  # Investigate top findings
  specs <- unique(unlist(sapply(top.findings, rownames)))[1:topN] # phylotypes

  if (length(specs)>0) {
    specs.shortnames <- sapply(specs, function(x) {ss <- strsplit(x, " "); if (length(ss[[1]]) > 1) { paste(paste(substr(ss[[1]][[1]], 1, 1), ".", sep = ""), paste(ss[[1]][-1], collapse = " "), sep = " ")} else {ss[[1]][[1]]}})
    names(specs.shortnames) <- specs
    s <- rownames(smat) # samples
    df <- cbind(smat[s, specs], annot[s,])
    means <- reshape::melt(aggregate(df[specs], by=list(varname = df[[varname]]), FUN=mean))
    stds <- reshape::melt(aggregate(df[specs], by=list(varname = df[[varname]]), FUN=sd))
    Nsqrt <- sqrt(as.numeric(reshape::melt(aggregate(df[specs], by = list(varname = df$varname), FUN = length))$value))
    dfm <- cbind(means, mean = means$value, sd = stds$value, sd.of.mean = as.numeric(stds$value) / Nsqrt)
    dfm$shortnames <- as.factor(specs.shortnames[dfm$variable])

    # Create the barplot component
    p <- ggplot2::ggplot(dfm, aes(x = shortnames, y = value, fill = varname)) 
    dodge <- position_dodge(width = 0.9)
    p <- p + geom_bar(position="dodge") 
    p <- p + geom_errorbar(aes(x = shortnames, ymax = mean + 1.96*sd.of.mean, ymin = mean - 1.96*sd.of.mean), position = dodge, width=0.25)
    p <- p + scale_fill_grey() + theme_bw() + ylab("Signal") + xlab("") 
    p <- p + opts(axis.text.x=theme_text(angle=-20)) 

    } else {
      warning("No features available.")
      p <- NULL  
    }

}



#' Description: Visualization of top findings from pairwise.comparisons function
#'              
#' Arguments:
#'   @param x data matrix
#'   @param y annotaction factor
#'   @param oligo.map mapping between features
#'   @param color.level feature level to color
#'   @param bar.level feature level for bars
#'   @param top.findings output from pairwise.comparisons function
#'   @param top.findings.qvals output from pairwise.comparisons function
#'   @param qth q-value threshold
#'   @param qth.star q-value threshold for stars
#'   @param mode barplot / heatmap
#' Returns:
#'   @return List with top findings from pairwise comparisons and their q-values
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

PlotTopComparisons <- function (x, y, oligo.map, color.level, bar.level, top.findings, top.findings.qvals, qth = 0.05, qth.star = 0.01, mode = "barplot") {

  # x <- smat; color.level = "level 1"; bar.level = "level 2"; qth <- 0.05; mode = "heatmap"

  flevels <- as.character(unique(y))

  map <- oligo.map[[color.level]]
  names(map) <- as.character(oligo.map[[bar.level]])

  # Pick the most significant findings only
  nams <- unique(unlist(sapply(top.findings, function(tab){unlist(rownames(subset(tab, qvalue < qth)))})))

  x <- x[, nams]  

  fc <- aggregate(x, by = list(varname = y), mean) # Mean for each group

  # Compare each group to other group
  comparisons <- list()

  for (i in 1:(length(flevels)-1)) {
    for (j in (i+1):length(flevels)) {
      di <- flevels[[i]]
      dj <- flevels[[j]]
      nam <- paste(di, dj, sep = "-")
      comparisons[[nam]] <- as.numeric(subset(fc, varname == di)[-1]) - as.numeric(subset(fc, varname == dj)[-1])
    }
  }

  nams <- names(comparisons)
  df <- data.frame(comparisons)
  colnames(df) <- nams
  rownames(df) <- colnames(x)
  df$bar.level <- colnames(x)
  df$color.level <- factor(map[as.character(df$bar.level)]) # Check L2 for the phylotype species

  if (mode == "barplot") {

    specs.shortnames <- sapply(as.character(df$bar.level), function(x) {ss <- strsplit(x, " "); if (length(ss[[1]]) > 1) { paste(paste(substr(ss[[1]][[1]], 1, 1), ".", sep = ""), paste(ss[[1]][-1], collapse = " "), sep = " ")} else {ss[[1]][[1]]}})
    names(specs.shortnames) <- as.character(df$bar.level)
    df$shortnames <- factor(specs.shortnames[df$bar.level])
    dfm <- reshape::melt(df)

    p <- ggplot2::ggplot(dfm) 
    dodge <- position_dodge(width = 0.9) # define the width of the dodge
    p <- p + geom_bar(position = "dodge", width = 0.4, aes(x = shortnames, y = value, fill = color.level))
    p <- p + facet_grid(.~variable, scales = "free", space = "free")
    p <- p + coord_flip()
    p <- p + theme_bw() + ylab("") + xlab("Fold Change")
    p <- p + opts(axis.text.x = theme_text(size = 7))
    # p <- p + opts(legend.position = c(0.15,.15))

  } else if (mode == "heatmap") {

    ncomp <- ncol(df) - 2 # number of comparisons

    dfm <- reshape::melt(df[, 1:(ncomp + 1)], variable = "comparison") # reshape::melt matrix
    keep <- !is.na(dfm$value)
    dfm <- dfm[keep, ]
    dfm$comparison <- droplevels(dfm$comparison)
    #flevels <- levels(dfm$comparison)

    qmat <- NULL; nams <- c()
    if (length(flevels)>1) {
      for (i in 1:(length(flevels)-1)) {
        for (j in (i+1):(length(flevels))) {
  
          ai <- flevels[[i]]
          aj <- flevels[[j]]
          nam <- paste(ai, aj, sep = "-")

          qmat <- cbind(qmat, unlist(top.findings.qvals[[ai]][[aj]])[df[[bar.level]]])
          nams <- c(nams, nam)
        }
      }
    } else {

          nam <- flevels[[1]]
          qmat <- top.findings.qvals[[1]][[1]][df$bar.level]
	  qmat <- matrix(qmat, length(qmat))
	  rownames(qmat) <- df$bar.level
          nams <- c(nams, nam)

    }

    if (nrow(qmat) > 0) {

      qmat <- data.frame(qmat)
      colnames(qmat) <- nams
 
      qmat[["phylotype"]] <- rownames(qmat)
      dfm.qval <- reshape::melt(qmat)
      names(dfm.qval) <- "value"

      cex.xlab = 7; cex.ylab = 12

      p <- ggplot2::ggplot(dfm, aes(x = comparison, y = bar.level, fill = value))
      p <- p + geom_tile() 
      limit <- ceiling(max(abs(na.omit(dfm$value))))
      p <- p + scale_fill_gradientn("Fold Change", breaks=seq(from=limit, to=-limit, by=-0.2), colours=c("darkblue", "blue", "white", "red", "darkred"), limits=c(-limit,limit)) 
      p <- p + theme_bw() 
      p <- p + opts(axis.text.x=theme_text(angle = 0, size = cex.xlab)) 
      p <- p + opts(axis.text.y=theme_text(size = cex.ylab))

      # Merkkaa merkitsevat tahdella
      p <- p + geom_text(data=subset(dfm.qval, value < qth.star), aes(x = variable, y = value, label = "+"), col = "white", size = 3)
      p <- p + xlab("") + ylab("")
    }
  }

  p

}

#' Description: ggplot2::ggplot2 theme
#'              
#' Arguments:
#'   @param colour colour
#'   @param size size
#'   @param linetype linetype
#' Returns:
#'   @return set theme
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

theme_bottom_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot2::ggplot(...) + opts( panel.border=theme_bottom_border() ) + ...
  structure(
    function(x = 0, y = 0, width = 1, height = 1, ...) {
      polylineGrob(
        x=c(x, x+width), y=c(y,y), ..., default.units = "npc",
        gp=gpar(lwd=size, col=colour, lty=linetype),
      )
    },
    class = "theme",
    type = "box",
    call = match.call()
  )
}