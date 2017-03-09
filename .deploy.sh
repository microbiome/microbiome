#!/bin/bash
# See https://medium.com/@nthgergo/publishing-gh-pages-with-travis-ci-53a8270e87db

set -e # Exit with nonzero exit code if anything fails

SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"

# Save some useful information
#SHA=`git rev-parse --verify HEAD`
GH_REPO="@github.com/microbiome/microbiome.git"
FULL_REPO="https://$GH_TOKEN$GH_REPO"

#function doCompile {
  # run pkgdown, put results in 'docs' directory,i
  # and don't paste the results of the examples
  # then copy the whole thing to `out`
#  Rscript -e "pkgdown::build_site(path = 'docs', examples = TRUE)"
#  cp -R ../docs .
#  ls ../ > files3.txt  
  #  ./compile.sh
#}

# Pull requests and commits to other branches shouldn't try to deploy, just build to verify
if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
    echo "Skipping deploy; just doing a build."
    exit 0
fi

rm -rf public # ; || exit 0;
mkdir -p public

# config
git config --global user.email "nobody@nobody.org"
git config --global user.name "Travis CI"

# Unzip the newly created R package 
tar -zxvf *.tar.gz

# Deploy
cd public

git init

# Copy the vignettes from the newly generated package in here
cp ../microbiome/inst/doc/*.html .
ls ../ > files5.txt

# Run our compile script
# doCompile
# ls ../ > files3.txt  
#Rscript -e "pkgdown::build_site(path = 'docs', examples = TRUE)"
#Rscript -e "library(pkgdown)"
R CMD BATCH "1" # pkgdown::build_site(path = 'docs', examples = TRUE)
#cp -R ../docs .
ls ../ > files7.txt

# Add to git and deploy
git add *.html
#git add docs
git add files*.txt
git commit -a -m "Deploy to Github Pages"
git push --force --quiet $FULL_REPO $SOURCE_BRANCH:$TARGET_BRANCH # > /dev/null 2>&1

# --------------------------------

#for files in '../microbiome/vignettes/*.html'; do
#        cp $files .
#done



