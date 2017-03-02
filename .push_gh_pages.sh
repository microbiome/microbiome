#!/bin/bash

# Modified from:
# https://rmflight.github.io/posts/2014/11/travis_ci_gh_pages.html

# Completely overwrites gh-pages
# would be also possible to pull and modify instead
rm -rf out; || exit 0;
mkdir out;

#SOURCE_BRANCH="master"
#TARGET_BRANCH="gh-pages"

# Can add separate vignette build here later
#function doCompile {
#  ./compile.sh
#}
# Run our compile script -> Add this to a suitable place
# doCompile


# Pull requests and commits to other branches shouldn't try to deploy, just build to verify
#if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
#    echo "Skipping deploy; just doing a build."
#    doCompile
#    exit 0
#fi

# Save some useful information
#REPO=`git config remote.origin.url`
#SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
#SHA=`git rev-parse --verify HEAD`
GH_REPO="@github.com/microbiome/microbiome.git"
FULL_REPO="https://$GH_TOKEN$GH_REPO"

for files in '*.tar.gz'; do
        tar xfz $files
done

cd out
git init
git config user.name "microbiome-travis"
git config user.email "travis"
#cp ../microbiome/inst/doc/vignette.html index.html
#cp ../microbiome/vignettes/vignette.html index.html
# Add when index works
#for files in '../optiRum/inst/doc/*.html'; do
#        cp $files .
#done

touch index.html

git add .
git commit -m "Deployed to github pages"
git push --force --quiet $FULL_REPO master:gh-pages

