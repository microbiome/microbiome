#!/bin/bash

# Modified from:
# https://rmflight.github.io/posts/2014/11/travis_ci_gh_pages.html

# Completely overwrites gh-pages
# would be also possible to pull and modify instead
rm -rf out || exit 0;
mkdir out;

GH_REPO="@github.com/microbiome/microbiome.git"

FULL_REPO="https://$GH_TOKEN$GH_REPO"

for files in '*.tar.gz'; do
        tar xfz $files
done

cd out
git init
#git config user.name "microbiome-travis"
git config user.name "antagomir"
git config user.email "travis"
#cp ../microbiome/inst/doc/vignette.html index.html
cp ../microbiome/vignettes/vignette.html index.html

git add .
git commit -m "deployed to github pages"
git push --force --quiet $FULL_REPO master:gh-pages
