# Copyright (C) 2011-2013 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Default list of removed phylotypes and oligos
#'
#' Arguments:
#'  @param chip Chip name (HIT/MIT/PIT/Chick)Chip
#' Returns:
#'   @return List of removed oligos and phylotypes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

phylotype.rm.list <- function (chip) {

  rm.phylotypes <- list()

  if (chip == "HITChip") {
    
    rm.phylotypes[["oligos"]] <- c("UNI 515", "HIT 5658", "HIT 1503", "HIT 1505", "HIT 1506")
    rm.phylotypes[["species"]] <- c("Victivallis vadensis")
    rm.phylotypes[["L1"]] <- c("Lentisphaerae")
    rm.phylotypes[["L2"]] <- c("Victivallis")

  } else if (chip == "MITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["L1"]] <- c()
    rm.phylotypes[["L2"]] <- c()

  } else if (chip == "PITChip") {

    # Based on JZ mail 9/2012; LL

    rm.old.oligos <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.new.oligos <- c("PIT_1083", "PIT_1022", "PIT_1057", "PIT_1023", "PIT_1118", "PIT_1040", "PIT_1058", "PIT_1119", "PIT_122", "PIT_1221", "PIT_1322", "PIT_1367", "PIT_1489", "PIT_160", "PIT_1628", "PIT_1829", "PIT_1855", "PIT_1963", "PIT_1976", "PIT_1988", "PIT_2002", "PIT_2027", "PIT_2034", "PIT_2101", "PIT_2196", "PIT_2209", "PIT_2281", "PIT_2391", "PIT_2392", "PIT_2418", "PIT_2425", "PIT_2426", "PIT_2498", "PIT_2555", "PIT_2563", "PIT_2651", "PIT_2654", "PIT_2699", "PIT_2741", "PIT_2777", "PIT_2786", "PIT_2936", "PIT_35", "PIT_425", "PIT_427", "PIT_428", "PIT_429", "PIT_435", "PIT_481", "PIT_605", "PIT_7", "PIT_733", "PIT_734", "PIT_892")
    rm.phylotypes[["oligos"]] <- c(rm.old.oligos, rm.new.oligos)
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["L0"]] <- c("Nematoda", "Apicomplexa", "Euryarchaeota", "Ascomycota", "Parabasalidea", "Chordata")
    rm.phylotypes[["L1"]] <- c("Chromadorea", "Coccidia", "Methanobacteria", "Saccharomycetales", "Trichomonada", "Mammalia")
    rm.phylotypes[["L2"]] <- c("Ascaris suum et rel.", "Eimeria  et rel.", "Methanobrevibacter et rel.", "Saccharomyces et rel.", "Trichomonas et rel.", "Uncultured Mammalia", "Uncultured methanobacteria")

  } else if (chip == "ChickChip") {
    warning("No universal probes excluded from ChichChip yet!")
  }

  rm.phylotypes

}
#' Description: species-levels mappings
#'
#' Arguments:
#'   @param spec species
#'   @param phylogeny.info phylogeny.info
#'
#' Returns:
#'   @return species-levels mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

species2levels <- function (spec, phylogeny.info) {

   phylogeny.info <- polish.phylogeny.info(phylogeny.info)
	       
   omap <- phylogeny.info[match(as.character(spec), phylogeny.info$species), c("species", "L2", "L1")]
   omap[["species"]] <- factor(omap[["species"]])
   omap[["L1"]] <- factor(omap[["L1"]])
   omap[["L2"]] <- factor(omap[["L2"]])

   omap
}

#' Description: level2-level1 mappings
#'
#' Arguments:
#'   @param l2 level2 phylotypes
#'   @param phylogeny.info phylogeny.info
#'
#' Returns:
#'   @return mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

level2TOlevel1 <- function (l2, phylogeny.info) {

   phylogeny.info <- polish.phylogeny.info(phylogeny.info)

   omap <- phylogeny.info[match(as.character(l2), phylogeny.info[["L2"]]), c("L2", "L1")]
   omap[["L2"]] <- factor(omap[["L2"]])
   omap[["L1"]] <- factor(omap[["L1"]])

   omap
}


#' Description: L2-L0 mappings
#'
#' Arguments:
#'   @param l2 level2 phylotypes
#'   @param phylogeny.info phylogeny.info
#'
#' Returns:
#'   @return mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

level2TOlevel0 <- function (l2, phylogeny.info) {

   phylogeny.info <- polish.phylogeny.info(phylogeny.info)

   omap <- phylogeny.info[match(as.character(l2), phylogeny.info[["L2"]]), c("L2", "L0")]
   omap[["L2"]] <- factor(omap[["L2"]])
   omap[["L0"]] <- factor(omap[["L0"]])

   omap
}


#' Description: L1-L0 mappings
#'
#' Arguments:
#'   @param l1 level1 phylotypes
#'   @param phylogeny.info phylogeny.info
#'
#' Returns:
#'   @return mappings
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

level1TOlevel0 <- function (l1, phylogeny.info) {

   phylogeny.info <- polish.phylogeny.info(phylogeny.info)

   omap <- phylogeny.info[match(as.character(l1), phylogeny.info[["L1"]]), c("L1", "L0")]
   omap[["L1"]] <- factor(omap[["L1"]])
   omap[["L0"]] <- factor(omap[["L0"]])

   omap
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
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

list.color.scales <- function () {
  ## Different colour scales
  list('white/blue'=colorRampPalette(c("white","darkblue"),interpolate='linear')(100),
       'white/black'=colorRampPalette(c("white","black"),interpolate='linear')(100),
       'black/yellow/white'=colorRampPalette(c("black","yellow","white"),bias=0.5,interpolate='linear')(100))

}



#' Description: Format string vector to mysql query format
#' 
#' Arguments:
#'   @param s string vector
#'
#' Returns:
#'   @return mysql query version
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

mysql.format <- function (s) {
  paste("('",paste(as.character(s),collapse="','",sep="'"),"')",sep="")
}




#' Description: Calculate species summaries and possibly update d.oligo2
#'
#' Arguments:
#'   @param d.oligo2 d.oligo2
#'   @param bgc.method background correction method
#' Returns:
#'   @return Background-corrected data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

oligo.bg.correction <- function (d.oligo2, bgc.method) {

  if ( bgc.method == "6*sd bkg intensity" ){ bgth <- 6 }

  d.oligo2 <- threshold.data(d.oligo2, bgth)
  d.oligo2 <- apply(d.oligo2, c(1,2), function(x) max(0, x))
  
  d.oligo2

}



#' Description: Check number of matching phylotypes for each probe
#' 
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param level phylotype level
#'
#' Returns:
#'   @return number of matching phylotypes for each probe
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

n.phylotypes.per.oligo <- function (phylogeny.info, level) {
  sapply(split(phylogeny.info[, c("oligoID", level)], phylogeny.info$oligoID), function(x) length(unique(x[[level]])))
}

#' Description: Write matrix in tab file
#'
#' Arguments:
#'   @param dat data matrix
#'   @param filename output file
#'   @param verbose verbose
#' Returns:
#'   @return output file location
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

WriteMatrix <- function (dat, filename, verbose = FALSE) { 

  if (verbose) { message(paste("Writing output in ", filename)) }
  write.table(dat, file = filename, quote = FALSE, sep = "\t", row.names = FALSE)
  filename

}



#' Description: determine threshold for bg correction
#'
#' Arguments:
#'   @param dat data matrix (in approximately normal scale ie. logged)
#'
#' Returns:
#'   @return threshold value
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate.min.threshold <- function (dat) {

  #estimate min threshold 
  DD <- density(as.numeric(unlist(dat)))

  #find mode
  noise_mode <- DD$x[[which.max(DD$y)]] # DD$x[which(diff(DD$y)<0)[1]]

  #compute sd of noise
  noise_sd <- sd(dat[dat < noise_mode])

  #threshold
  low.thresh <- noise_mode + 6*noise_sd

  low.thresh
}


#' Description: determine detection threshold for the data
#'
#' Arguments:
#'   @param dat data
#'   @param sd.times standard deviation threshold
#'
#' Returns:
#'   @return thresholded data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

threshold.data <- function(dat, sd.times = 6){

  thr <- apply(dat, 2, function(x){
      DD <- density(as.numeric(x), adjust = 1.2, na.rm = T);
      noise_mode <- DD$x[which(DD$y==max(DD$y))[1]];
      noise_sd   <- sd(x[x < noise_mode], na.rm = T);
      low.thresh <- noise_mode + sd.times*noise_sd;
      low.thresh 
    })

  # Subtract background from signal intensities in each sample
  data.mat <- t(apply(dat, 1, function(Tr){ Tr-thr })) 
  return(data.mat)
}



#' Description: Probeset summarization with various methods.
#' 
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose print intermediate messages
#'   @param rm.species Species to exclude
#' Returns:
#'   @return summarized data matrix in log10 scale
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.species <- function (phylogeny.info, oligo.data, method, verbose = TRUE, rm.species = c("Victivallis vadensis")) {

  level <- "species"			    

  if (method == "nmf") {warning("nmf oligo summarization method not implemented at species level"); return(NULL)}
  probesets <- retrieve.probesets(phylogeny.info, level = level)
  probesets <- probesets[setdiff(names(probesets), rm.species)]
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets), ncol(oligo.data)), 
  		    	      dimnames = list(names(probesets), colnames(oligo.data))) 

  for (set in names(probesets)) {

    if (verbose) { message(set) }

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "rpa") {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      vec <- rpa.fit(dat, tau2.method = "robust", alpha = 1 + 0.1*ncol(oligo.data)/2, beta = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2)$mu

    } else if (method == "ave") {

      vec <- log10(colMeans((10^dat), na.rm = T))

    } else if (method == "sum") {

      # Weight each probe by the inverse of the number of matching phylotypes
      # Then calculate sum -> less specific probes are downweighted
      # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
      dat2 <- (10^dat) / nPhylotypesPerOligo[rownames(dat)]
      dat2[dat2 < 1] <- 1
      vec <- log10(colSums(dat2, na.rm = T))

    }
    
    summarized.matrix[set, ] <- vec 

  }

  summarized.matrix
  
}

#' Description: RPA for HITChip
#' 
#' Arguments:
#'   @param level level
#'   @param phylo phylo
#'   @param oligo.data oligo.data
#'
#' Returns:
#'   @return RPA preprocessed data
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

calculate.rpa <- function (level, phylo, oligo.data) {

  # List entities (e.g. species)
  phylo.list <- split(phylo, phylo[[level]])
  entities <- names(phylo.list)

  # initialize
  summarized.matrix <- array(NA, dim = c(length(entities), ncol(oligo.data)), dimnames = list(entities, colnames(oligo.data)))
  noise.list <- list() 

  for (entity in names(phylo.list)) {
    message(entity)

    # Pick expression for particular probes
    probes <- unique(phylo.list[[entity]][, "oligoID"])

    # oligo.data is already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 

    # dat: probes x samples
    if (nrow(dat) < 2) {
      vec <- as.vector(dat) # NOTE: circumvent RPA if there are no replicates 
      noise <- NA
    } else {
      res <- rpa.fit(dat)
      vec <- res$mu
      noise <- sqrt(res$tau2)
      names(noise) <- probes
    }

    noise.list[[entity]] <- noise
    summarized.matrix[entity, ] <- vec #, epsilon, alpha, beta, tau2.method, d.method)

  }

  list(emat = summarized.matrix, noise = noise.list)
  
}



