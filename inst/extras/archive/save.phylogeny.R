message("Fetching Phylogeny from the database")
library(microbiome)
library(HITChipDB)
dbuser = "root"; dbpwd = "fidipro"; dbname = "phyloarray"; host = '127.0.0.1'; port = 3307

message("Fetching Phylogeny from the database")
phylogeny.full <- get.phylogeny.info("16S", 
	    		     dbuser = dbuser, 
			     dbpwd = dbpwd, 
			     dbname = dbname, 
			     host = host, 
			     port = port,
			     chip = "HITChip")

phylogeny.filtered.orig <- GetPhylogeny("HITChip", "filtered")

# Fix taxon name
phylogeny.full$L2 <- gsub("^Clostridia$", "Clostridium (sensu stricto)", as.character(phylogeny.full$L2))

# This handles also pmTm, complement and mismatch filtering
# This is the phylogeny used in probe summarization into taxonomic levels
rm.phylotypes <- phylotype.rm.list("HITChip") 
rm.oligos <- sync.rm.phylotypes(rm.phylotypes, phylogeny.full)$oligos

phylogeny.filtered <- prune16S(phylogeny.full, pmTm.margin = 2.5, complement = 1, mismatch = 0, rmoligos = rm.oligos, remove.nonspecific.oligos = FALSE) # remove.nonspecific.oligos refers to probes with multiple L2 groups

# Remove probes that target multiple L1 groups
lmap <- levelmap(NULL, level.from = "oligoID", level.to = "L1", phylogeny.info = phylogeny.filtered)
hits <- sapply(lmap, length)
ambiguous.l1.probes <- unique(names(which(hits > 1)))

phylogeny.filtered <- subset(phylogeny.filtered, !oligoID %in% ambiguous.l1.probes)

# Keep only relevant cols
phylogeny.full <- phylogeny.full[, 1:6]; 
phylogeny.filtered <- phylogeny.filtered[, 1:6]; 

# Remove duplicate rows
phylogeny.full <- phylogeny.full[!duplicated(phylogeny.full),]
phylogeny.filtered <- phylogeny.filtered[!duplicated(phylogeny.filtered),]
 
# Write to file
write.table(phylogeny.filtered, file = "../extdata/phylogeny.filtered.tab", quote = F, row.names = F, sep = "\t")

