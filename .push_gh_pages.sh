#!/bin/bash

rm -rf out || exit 0;
mkdir out;

GH_REPO="@github.com/microbiome/microbiome.git"

FULL_REPO="https://$GH_TOKEN$GH_REPO"

for files in '*.tar.gz'; do
        tar xfz $files
done

cd out
git init
git config user.name "microbiome-travis"
git config user.email "leo.lahti@iki.fi"
cp ../microbiome/inst/doc/vignette.html index.html

git add .
git commit -m "deployed to github pages"
git push --force --quiet $FULL_REPO master:gh-pages
