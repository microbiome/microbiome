# Mahdollisesti riittaa:
# sudo apt-get update
# sudo apt-get install libxml2-dev
# sudo apt-get install emacs24
# sudo apt-get install git
# sudo apt-get upgrade r-base
# cd ~/bin/
# wget http://cran-mirror.cs.uu.nl/src/base/R-3/R-3.2.0.tar.gz
# tar -zxvf R-3.2.0.tar.gz
# ./configure
# make
# sudo make install

# In R:
install.packages("devtools")
source("http://www.bioconductor.org/biocLite.R")
biocLite("WGCNA")
biocLite("BiocStyle")
biocLite("knitr")
biocLite("rmarkdown")
biocLite("phyloseq")
biocLite("RPA")

library(devtools)
install_github("microbiome/microbiome")
install_github("microbiome/HITChipDB")

git config --global user.email "leo.lahti@iki.fi"
git config --global user.name "antagomir"

