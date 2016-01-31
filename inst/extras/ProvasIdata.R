# This script retrieves the data from 
# Lahti et al. PeerJ 1:e32, 2013 https://peerj.com/articles/32/
# and saves it into a more compact matrix format.

# Download and extract the tarball (data stored in Data/ directory)::
system("wget https://dfzljdn9uc3pi.cloudfront.net/2013/32/1/SupplementalDataS2.tar.gz")
system("tar -zxvf SupplementalDataS2.tar.gz")

# Read HITChip genus-level intestinal microbiota-profiling data
hit <- read.csv("Data/SupplementalDataHITChipS3.tab", sep = " ", skip = 1, header = FALSE, row.names = 1)

# The taxon names are unfortunately space-separated, so they need to be
# specified and read in manually from the file:
h <- unlist(strsplit(readLines("Data/SupplementalDataHITChipS3.tab", n = 1), " "))
inds <- list(1,2,3,4,5:8,9:11,12,13,14:17,18:21,22:25,26,27,28:30,31,32,33:36,37:40,41:44,45:48,49:52,53:56,57:60,61:64,65,66:68,69,70:73,74:77,78,79:82,83,84:87,88,89:92,93:96,97:100,101:104,105:108,109:112,113:116,117:120,121:124,125:128,129:132,133:136,137,138:141,142:145,146,147:149,150,151:154,155:158,159:162,163,164:167,168:171,172:175,176:179,180:183,184:187,188:191,192:195,196:199,200,201,202,203,204,205:208,209:212,213:216,217:220,221:224,225:228,229:232,233,234,235:238,239:242,243,244,245:248,249,250,251,252:255,256:259,260:263,264:267,268:271,272:275,276:279,280:283,284:287,288:291,292:295,296:299,300:303,304,305:307,308,309:312,313:316,317:320,321:324,325:328,329:332,333,334:337,338,339:342,343:346,347:350,351:354,355:358,359:361,362:363,364:365,366:368,369:371,372:373,374:375,376,377:379,380,381:383,384,385:387)
nams <- sapply(inds, function (i) {paste(h[i], collapse = " ")})

# Add the retrieved taxon names to HITChip data columns
colnames(hit) <- nams

# Read Lipidomic profiling data
lip <- read.csv("Data/SupplementalDataLipidsS2.tab", sep = " ", row.names = 1)
lip2 <- apply(lip, 2, function (x) as.numeric(gsub("\\,", "\\.", x)))
dimnames(lip2) <- dimnames(lip)

# Write to file
write.table(hit, file = "~/Rpackages/microbiome/microbiome/inst/extdata/ProvasI-L2.tab", sep = ";", quote = FALSE)
write.table(lip2, file = "~/Rpackages/microbiome/microbiome/inst/extdata/ProvasI-Lipid.tab", sep = ";", quote = FALSE)

# ........................

data(peerj32)
# lipids microbes meta
# Taxonomy
taxonomy <- GetPhylogeny("HITChip", "filtered")
taxonomy <- unique(as.data.frame(taxonomy[, c("L1", "L2")]))
rownames(taxonomy) <- as.vector(taxonomy[, "L2"])

# Merging data matrices into phyloseq format:
pseq <- hitchip2physeq(t(otu), meta, taxonomy)

round(peerj32$microbes - min(peerj32$microbes))