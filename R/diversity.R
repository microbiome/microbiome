#' estimate.diversity
#'
#' Description: Estimate diversities for each sample (column)
#' Aliases: get.diversity.estimates. Also replaces the function ST: diversity.indices
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param diversity.index diversity index
#'   @param det.th detection threshold. Used for richness and evenness estimation. Not used in diversity estimation. 
#'
#' Returns:
#'   @return diversity indicators
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

estimate.diversity <- function (dat, diversity.index = "shannon", det.th = NULL) {

  veganT <- require(vegan)
  if(!veganT) { install.packages("vegan") }

  # Specify detection threshold if not provided
  if (is.null(det.th)) {
    dat <- impute(dat) # impute the very few (3) missing values
    det.th <- 10^estimate.min.threshold(log10(dat))
    warning(paste("Applying detection threshold: ", det.th))
  }

  # Apply detection threshold
  dat <- dat - det.th
  dat[dat < 0] <- 0

  # Species diversity - only for species that exceeded the thresholding above
  H <- diversity(dat, index = diversity.index, MARGIN = 2)
  
  # Species richness - count phylotypes that exceed detection threshold
  S <- colSums(dat > 0)
  
  # Pielou's evenness (J) 
  H.shannon <- diversity(dat, index = "shannon", MARGIN = 2)
  J <- H.shannon/log(S)
  
  names(J) <- names(S) <- names(H) <- colnames(dat)

  data.frame(list(evenness = J, richness = S, diversity = H, det.th = det.th))
	
}

  # Compare to ACE estimate with different discretization bins
  # -> no clear way to discretize, will affect the results quite much
  # -> skip so far
  #par(mfrow=c(3,3))
  #for (discretization.resolution in 10^(-seq(-2,6,length=9))) {
  #  ab <- make.abundancy.table(dat, det.th, discretization.resolution)
  #  S2 <- estimateR(ab)["S.ACE", ]
  #  plot(S, S2, main = discretization.resolution)
  #}
  # Chao and ACE estimators by modifying this: requiring counts, add later?
  # Estimater(floor(10^t(dat))) 
  #S <- unlist(mclapply(1:ncol(dat), function(k) {'  
  #ab <- make.abundancy.table(dat[, k], discretization.resolution = 1e-4)
    #ChaoLee1992(ab, t=10, method="ACE")$Nhat
  #}))



#' make.abundancy.table
#'
#' Description: Calculate abundancies
#' Discretize Hitchip matrix to form abundancy table
#' of form j, nj where j is number of counts and nj is number
#' of phylotypes with the corresponding counts
#' this format is often required by richness estimation
#'
#' Arguments:
#'   @param dat data matrix
#'   @param det.th detection threshold
#'   @param discretization.resolution discretization resolution
#'
#' Returns:
#'   @return abundancy table
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

make.abundancy.table <- function (dat, det.th, discretization.resolution = 1) {

  di <- 10^dat - 10^det.th
  di[di<0] <- 0
  di <- discretization.resolution*(di)
  di[di>0 & di<1] <- 1
  di <- round(di)

  ab <- t(di)

  ab
}

#' PlotDiversity
#'
#' Description: Plot and save diversity indices
#' Based on ST 01-03-2011
#'
#' Arguments:
#'   @param dat data matrix (phylotypes x samples) in original (non-log) scale
#'   @param diversity.index Diversity index: "shannon" / "invsimpson"
#'   @param filename output file name to store the diversity indices
#'
#' Returns:
#'   @return diversity
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

PlotDiversity <- function (dat, diversity.index, filename = NULL) {

  div <- estimate.diversity(dat, diversity.index)$diversity

  ## Grouping:
  # Select samples
  samples <- colnames(dat)

  # Set the amount of groups you want to make
  groups <- c(1:4)	# Change the '4' into a higher number if you need/want more groups!
  G <- tk_select.list(groups, multiple = F, title = "Number of groups over which to divide your samples?")

  # Fill the groups with samples and give names
  ###### name all groups through user input
  ###### and allows the user to select the boxplot colours
  temp.color.new <- "" # needed later...
  for (i in 1:G) {
	# Group selecting
	title.g <- paste("Select samples belonging to group ",i,sep="")
	Gtemp <- tk_select.list(samples, multiple=T, 
	title=title.g)

	# Group naming (Not basic stuff here!)
	tt <- tktoplevel()
	name.g <- paste("Type the name of group ",i,sep="")
	tkwm.title(tt,name.g)
	textEntryVarTcl <- tclVar("group x")
	textEntryWidget <- tkentry(tt,width=40,textvariable=textEntryVarTcl)
	tkgrid(tklabel(tt,text="       "))
	tkgrid(tklabel(tt,text="Name:   "),textEntryWidget)
	tkgrid(tklabel(tt,text="       "))
	onOK <- function()
		{
		tempname <- tclvalue(textEntryVarTcl)
		tkdestroy(tt)
		}
	OK.but <- tkbutton(tt,text="   OK   ",command=onOK)
	tkgrid(OK.but)
	tkwait.window(tt)
	if(i==1){gr.list<-list(tempname = Gtemp)} else { gr.list <- c(gr.list,list(tempname=Gtemp)) }
	eval(parse(text=paste("names(gr.list)[i]<-tempname")))

	# Group colour selection (Not basic stuff here!)
	tc <- tktoplevel()
	tkwm.title(tc,paste("Color Selection for ",tempname,sep=""))
	color <- "blue"
	temp.color.new <<- color
	canvas <- tkcanvas(tc,width="80",height="25",bg=color)
	ChangeColor <- function() {
   		color <- tclvalue(tcl("tk_chooseColor",initialcolor=color,title="Choose a color"))
  		if (nchar(color)>0)
    		tkconfigure(canvas,bg=color)
		temp.color.new <<- color
		#cat("\n\nColour: ",color,"\n\n"); flush.console()
		}
	ChangeColor.button <- tkbutton(tc,text="Change Color",command=ChangeColor)
	onOK2 <- function() {
		tkdestroy(tc)
		}
	Ok.button <- tkbutton(tc,text="Ok, use this colour!",command=onOK2)
	tkgrid(canvas,ChangeColor.button,Ok.button)
	tkwait.window(tc)
	if(i==1){col.list<-list(tempname=temp.color.new)} else { col.list<-c(col.list,list(tempname=temp.color.new)) }
	eval(parse(text=paste("names(col.list)[i]<-tempname")))
  }

  ## "Box-plotting":
  if (diversity.index == "invsimpson") {diversity.text = "Inverse Simpson"}  
  if (diversity.index == "shannon") {diversity.text = "Shannon Diversity"}  
  boxplot(div[gr.list[[1]]], ylab = diversity.text, col=col.list[[1]],
  ylim=c(min(div),max(div)),xlim=c(0,length(gr.list))+0.5)
  axis(1,labels=names(gr.list),at=c(1:length(gr.list)))
  for (i in 2:length(gr.list)) {
    boxplot(div[gr.list[[i]]],col=col.list[[i]],names=names(gr.list)[i],add=T,at=i)
  }


  ## Step 3: Save tables
  save <- FALSE
  if (is.null(filename)) {
    save <- tk_select.list(c("yes","no"),multiple=F, title= "Save diversity data?")
    if(save == "yes") {
      filename <- choose.files(multi = F)
      save <- TRUE
    } 
  }

  if(save) {
    write.table(div, file = filename, sep="\t")
  }

  div

}

