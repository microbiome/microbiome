# https://support.rstudio.com/hc/en-us/articles/200626518-Customizing-Package-Build-Options
#
# Build, check and install the package from a clean tarball.
# Run this from inst/extras. Regenerate the man pages and README first:
#
#   $R/bin/Rscript document.R
#
# The version is read from DESCRIPTION so it does not need updating here.

# Override to check against another R, e.g. the one Bioc devel builds on:
#   R=$HOME/bin/R-devel/bin/R sh build.sh
R=${R:-$HOME/bin/R-4.5.1/bin/R}
TARBALL=microbiome_$(awk '/^Version:/ {print $2}' ../../DESCRIPTION).tar.gz

$R CMD build ../../ || exit 1 #--resave-data #--no-examples  --no-build-vignettes
$R CMD check $TARBALL #--no-build-vignettes --no-examples
$R --vanilla -q -e "BiocCheck::BiocCheck('$TARBALL')" # R CMD BiocCheck is gone
$R CMD INSTALL $TARBALL
