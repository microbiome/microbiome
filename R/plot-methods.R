#' @title Ordination Plot
#' @description General visualization for ordination.
#' @details Modified from phyloseq::plot_ordination to add some features.
#'          For examples of phyloseq ordination graphics, see
#' \href{http://joey711.github.io/phyloseq/plot_ordination-examples}{phyloseq online tutorials}.
#' @param physeq (Required). \code{\link{phyloseq-class}}. 
#'  The data about which you want to plot and annotate the ordination.
#' @param ordination (Required). An ordination object. Many different classes
#'  of ordination are defined by \code{R} packages. Ordination classes
#'  currently supported/created by the \code{\link{ordinate}} function are
#'  supported here. There is no default, as the expectation is that the 
#'  ordination will be performed and saved prior to calling this plot function.
#' @param type (Optional). The plot type. Default is \code{"samples"}. The
#'  currently supported options are
#'  \code{c("samples", "sites", "species", "taxa", "biplot", "split", "scree")}.
#'  The option ``taxa'' is equivalent to ``species'' in this case, and
#'  similarly, ``samples'' is equivalent to ``sites''.  The options
#'  \code{"sites"} and \code{"species"} result in a single-plot of just the 
#'  sites/samples or species/taxa of the ordination, respectively.
#'  The \code{"biplot"} and \code{"split"} options result in a combined
#'  plot with both taxa and samples, either combined into one plot (``biplot'')
#'  or separated in two facet panels (``split''), respectively.
#'  The \code{"scree"} option results in a call to \code{\link{plot_scree}},
#'  which produces an ordered bar plot of the normalized eigenvalues
#'  associated with each ordination axis. 
#' @param axes (Optional). A 2-element vector indicating the axes of the 
#'  ordination that should be used for plotting. 
#'  Can be \code{\link{character-class}} or \code{\link{integer-class}},
#'  naming the index name or index of the desired axis for the horizontal 
#'  and vertical axes, respectively, in that order. The default value, 
#'  \code{c(1, 2)}, specifies the first two axes of the provided ordination.
#' @param color (Optional). Default \code{NULL}. Character string.
#'  The name of the variable to map to
#'  colors in the plot. 
#'  This can be a sample variable 
#'  (among the set returned by \code{sample_variables(physeq)} )
#'  or taxonomic rank (among the set returned by \code{rank_names(physeq)}).
#'  The color scheme is chosen automatically
#'  by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_color_manual}}.
#' @param shape (Optional). Default \code{NULL}. Character string.
#'  The name of the variable to map
#'  to different shapes on the plot. 
#'  Similar to \code{color} option, but for the shape if points.
#'  The shape scale is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_shape_manual}}.
#' @param label (Optional). Default \code{NULL}. Character string.
#'  The name of the variable to map to text labels on the plot.
#'  Similar to \code{color} option, but for plotting text.
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#' @param show.density Show point density
#' @return A \code{\link{ggplot}} plot object, graphically summarizing
#'  the ordination result for the specified axes.
#' @seealso Many related examples are included in the phyloseq plot_ordination
#' page \href{http://joey711.github.io/phyloseq/plot_ordination-examples}{phyloseq online tutorials}.
#' @export
#' @examples 
#' # See other examples at
#' # http://joey711.github.io/phyloseq/plot_ordination-examples
#'   data(GlobalPatterns)
#'   taxa <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50])
#'   GP <- prune_taxa(taxa, GlobalPatterns)
#'   gp_bray_pcoa <- ordinate(GP, "CCA", "bray")
#'   plot_ordn(GP, gp_bray_pcoa, "samples", color="SampleType")
plot_ordn <- function(physeq, ordination, type="samples", axes=1:2,
                      color=NULL, shape=NULL, label=NULL, title=NULL,
		      show.density = TRUE){

  x <- y <- size <- ..density.. <- update_labels <- extract_eigenvalue <- rm.na.phyloseq <- NULL

  # Check input validity			    
  check <- plot_ordn_inputcheck(physeq, type, color, shape, label)
  type <- check$type
  color <- check$color
  shape <- check$shape
  label <- check$label

  if( type %in% c("scree") ){
    # Stop early by passing to plot_scree() if "scree" was chosen as a type
    return( phyloseq::plot_scree(ordination, title=title) )
  }

  # Prepare the data.frame
  DF <- plot_ordn_prepare_df(ordination, axes, physeq, type, color, shape)

  tmp <- check_variable_availability(DF, color, shape, label)
  color <- tmp$color
  shape <- tmp$shape
  label <- tmp$label

  # Prepare ordination map aes
  ord_map <- ordination_map(DF, type, color, shape, size) 

  # Construct the Plot 
  p <- ggplot(DF, ord_map)

  # Point density
  if (show.density) {
    bw.adjust <- 1
    bw <- bw.adjust * c(bandwidth.nrd(DF[,1]), bandwidth.nrd(DF[,2]))
    p <- p + stat_density2d(aes(fill = ..density.., color = NULL, shape = NULL), geom = "raster", h = 1, contour = FALSE)
    p <- p + scale_fill_gradient(low = "white", high = "black")
  }

  # Add points
  p <- p  + geom_point(na.rm=TRUE)

  # split/facet color and shape can be anything in one or other.
  if (type == "split") {
    # split-option requires a facet_wrap
    p <- p + facet_wrap(~id.type, nrow=1)
  }
  
  # If biplot, adjust scales
  if( type=="biplot" ){	
    if( is.null(color) ){
      # Rename color title in legend.
      p <- update_labels(p, list(colour="Ordination Type")) 
    } 
    # Adjust size so that samples are bigger than taxa by default.
    p <- p + scale_size_manual("type", values=c(Samples=5, Taxa=2))
  }
  
  # Add text labels to points
  if( !is.null(label) ){
    label_map <- aes_string(x=x, y=y, label=label, na.rm=TRUE)
    p = p + geom_text(label_map, data=rm.na.phyloseq(DF, label),
                      size=2, vjust=1.5, na.rm=TRUE)
  }
  
  # Add title
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }

  # FIXME extract_eigenvalue is missing often and this fails then
  # incuding the example case.
  skip <- TRUE # extract_eigenvalue is missing
  if (!skip) {
  # Add fraction variability to axis labels, if available
  if(length(extract_eigenvalue(ordination)[axes]) > 0 ){
    # Only attempt to add fraction variability
    # if extract_eigenvalue returns something
    eigvec <- extract_eigenvalue(ordination)
    # Fraction variability, fracvar
    fracvar = eigvec[axes] / sum(eigvec)
    # Percent variability, percvar
    percvar = round(100*fracvar, 1)
    # The string to add to each axis label, strivar
    # Start with the curent axis labels in the plot
    strivar = as(c(p$label$x, p$label$y), "character")
    # paste the percent variability string at the end
    strivar = paste0(strivar, "   [", percvar, "%]")
    # Update the x-label and y-label
    p <- p + xlab(strivar[1]) + ylab(strivar[2])
  }
  }
  
  # Return the ggplot object
  return(p)

}



check_variable_availability <- function (DF, color, shape, label) {

  # Check variable availability 
  if(!is.null(color)){ 
    if(!color %in% names(DF)){
      warning("Color variable was not in the data you provided.",
              "No color mapped.")
      color <- NULL
    }
  }
  if(!is.null(shape)){ 
    if(!shape %in% names(DF)){
      warning("Shape variable not in the data.",
              "No shape mapped.")
      shape <- NULL
    }
  }
  if(!is.null(label)){ 
    if(!label %in% names(DF)){
      warning("Label variable not in the data.",
              "No label mapped.")
      label <- NULL
    }
  }
  list(color = color, shape = shape, label = label)
}



ordination_map <- function (DF, type, color, shape, size) {

  # Grab the ordination axis names from the plot data frame (as strings)
  x = colnames(DF)[1]
  y = colnames(DF)[2]
  
  # Mapping section
  if( ncol(DF) <= 2){
    # If there is nothing to map, enforce simple mapping.
    message("No available covariate data to map on the points for this plot `type`")
    ord_map = aes_string(x=x, y=y)
  } else if( type %in% c("sites", "species", "split") ){
    ord_map = aes_string(x=x, y=y, color=color, shape=shape, na.rm=TRUE)
  } else if(type=="biplot"){
    # biplot, `id.type` should try to map to color and size.
    # Only size if color specified.
    if( is.null(color) ){
      ord_map = aes_string(x=x, y=y, size="id.type", color="id.type", shape=shape, na.rm=TRUE)
    } else {
      ord_map = aes_string(x=x, y=y, size="id.type", color=color, shape=shape, na.rm=TRUE)
    }
  }

  ord_map
  
}

plot_ordn_inputcheck <- function (physeq, type, color, shape, label) {

  if(length(type) > 1){
    warning("`type` can only be a single option,
            but more than one provided. Using only the first.")
    type <- type[[1]]
  }
  if(length(color) > 1){
    warning("The `color` variable argument should have length equal to 1.",
            "Taking first value.")
    color = color[[1]][1]
  }
  if(length(shape) > 1){
    warning("The `shape` variable argument should have length equal to 1.",
            "Taking first value.")
    shape = shape[[1]][1]
  }
  if(length(label) > 1){
    warning("The `label` variable argument should have length equal to 1.",
            "Taking first value.")
    label = label[[1]][1]
  }
  official_types = c("sites", "species", "biplot", "split", "scree")
  if(!inherits(physeq, "phyloseq")){
    if(inherits(physeq, "character")){
      if(physeq=="list"){
        return(official_types)
      }
    } 
    warning("Full functionality requires `physeq` be phyloseq-class ",
            "with multiple components.")
  }
  # Catch typos and synonyms
  type = gsub("^.*site[s]*.*$", "sites", type, ignore.case=TRUE)
  type = gsub("^.*sample[s]*.*$", "sites", type, ignore.case=TRUE)
  type = gsub("^.*species.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*taxa.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*OTU[s]*.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*biplot[s]*.*$", "biplot", type, ignore.case=TRUE)
  type = gsub("^.*split[s]*.*$", "split", type, ignore.case=TRUE)
  type = gsub("^.*scree[s]*.*$", "scree", type, ignore.case=TRUE)
  # If type argument is not supported...
  if( !type %in% official_types ){
    warning("type argument not supported. `type` set to 'samples'.\n",
            "See `plot_ordination('list')`")
    type <- "sites"
  }

  list(type = type, color = color, shape = shape, label = label)

}




plot_ordn_prepare_df <- function (ordination, axes, physeq, type, color, shape) {

  veganifyOTU <- access <- rp.joint.fill <- NULL

  # The plotting data frames.
  # Call scores to get coordinates.
  # Silently returns only the coordinate systems available.
  # e.g. sites-only, even if species requested.
  specDF = siteDF = NULL
  trash1 = try({siteDF <- scores(ordination, choices = axes, 
                                 display="sites", physeq=physeq)},
               silent = TRUE)
  trash2 = try({specDF <- scores(ordination, choices = axes, 
                                 display="species", physeq=physeq)},
               silent = TRUE)
  # Check that have assigned coordinates to the correct object
  siteSampIntx = length(intersect(rownames(siteDF), sample_names(physeq)))
  siteTaxaIntx = length(intersect(rownames(siteDF), taxa(physeq)))
  specSampIntx = length(intersect(rownames(specDF), sample_names(physeq)))
  specTaxaIntx = length(intersect(rownames(specDF), taxa(physeq)))
  if(siteSampIntx < specSampIntx & specTaxaIntx < siteTaxaIntx){
    # Double-swap
    co = specDF
    specDF <- siteDF
    siteDF <- co
    rm(co)
  } else {
    if(siteSampIntx < specSampIntx){
      # Single swap
      siteDF <- specDF
      specDF <- NULL
    }
    if(specTaxaIntx < siteTaxaIntx){
      # Single swap 
      specDF <- siteDF
      siteDF <- NULL
    }
  }

  # Define a function to check if a data.frame is empty
  is_empty = function(x){
    length(x) < 2 | suppressWarnings(all(is.na(x)))
  }
  
  # If both empty, warn and return NULL
  if(is_empty(siteDF) & is_empty(specDF)){
    warning("Could not obtain coordinates from the provided `ordination`. \n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  # If either is missing, do weighted average
  if(is_empty(specDF) & type != "sites"){
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    specDF <- data.frame(wascores(siteDF, w = veganifyOTU(physeq)), stringsAsFactors=FALSE)
  }
  if(is_empty(siteDF) & type != "species"){ 
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    siteDF <- data.frame(wascores(specDF, w = t(veganifyOTU(physeq))), stringsAsFactors=FALSE)
  }
  # Double-check that have assigned coordinates to the correct object
  specTaxaIntx <- siteSampIntx <- NULL
  siteSampIntx <- length(intersect(rownames(siteDF), sample_names(physeq)))
  specTaxaIntx <- length(intersect(rownames(specDF), taxa(physeq)))
  if(siteSampIntx < 1L & !is_empty(siteDF)){
    # If siteDF is not empty, but it doesn't intersect the sample_names in physeq, warn and set to NULL
    warning("`Ordination site/sample coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    siteDF <- NULL
  }
  if(specTaxaIntx < 1L & !is_empty(specDF)){
    # If specDF is not empty, but it doesn't intersect the taxa in physeq, warn and set to NULL
    warning("`Ordination species/OTU/taxa coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    specDF <- NULL
  }
  # If you made it this far and both NULL, return NULL and throw a warning
  if(is_empty(siteDF) & is_empty(specDF)){
    warning("Could not obtain coordinates from the provided `ordination`. \n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  if(type %in% c("biplot", "split") & (is_empty(siteDF) | is_empty(specDF)) ){
    # biplot and split require both coordinates systems available. 
    # Both were attempted, or even evaluated by weighted average.
    # If still empty, warn and switch to relevant type.
    if(is_empty(siteDF)){
      warning("Could not access/evaluate site/sample coordinates. Switching type to 'species'")
      type <- "species"
    }
    if(is_empty(specDF)){
      warning("Could not access/evaluate species/taxa/OTU coordinates. Switching type to 'sites'")
      type <- "sites"
    }
  }
  if(type != "species"){
    # samples covariate data frame, `sdf`
    sdf <- NULL
    sdf <- data.frame(sample_data(physeq), stringsAsFactors=FALSE)
    if( !is_empty(sdf) & !is_empty(siteDF) ){
      # The first two axes should always be x and y, the ordination axes.
      siteDF <- cbind(siteDF, sdf[rownames(siteDF), ])
    }
  }
  if(type != "sites"){
    # taxonomy data frame `tdf`
    tdf = NULL
    tdf = data.frame(tax_table(physeq), stringsAsFactors=FALSE)
    if( !is_empty(tdf) & !is_empty(specDF) ){
      # The first two axes should always be x and y, the ordination axes.
      specDF = cbind(specDF, tdf[rownames(specDF), ])
    }
  }
  # In "naked" OTU-table cases, `siteDF` or `specDF` could be matrix.
  if(!inherits(siteDF, "data.frame")){
    siteDF <- as.data.frame(siteDF, stringsAsFactors = FALSE)
  }  
  if(!inherits(specDF, "data.frame")){
    specDF <- as.data.frame(specDF, stringsAsFactors = FALSE)
  }


  # Define the main plot data frame, `DF`
  DF = NULL
  DF <- switch(EXPR = type, sites = siteDF, species = specDF, {
    # Anything else. In practice, type should be "biplot" or "split" here.
    # Add id.type label
    specDF$id.type <- "Taxa"
    siteDF$id.type <- "Samples"
    # But what if the axis variables differ b/w them?
    # Coerce specDF to match samples (siteDF) axis names
    colnames(specDF)[1:2] <- colnames(siteDF)[1:2]
    # Merge the two data frames together for joint plotting.
    DF = merge(specDF, siteDF, all=TRUE)
    # Replace NA with "samples" or "taxa", where appropriate (factor/character)
    if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "Samples") }
    if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "Taxa") }
    if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "Samples") }
    if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "Taxa") }
    DF
  })


  DF
}

