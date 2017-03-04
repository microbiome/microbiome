#!/bin/bash
# See https://medium.com/@nthgergo/publishing-gh-pages-with-travis-ci-53a8270e87db

#set -o errexit # put back

rm -rf public
mkdir public

# config
git config --global user.email "nobody@nobody.org"
git config --global user.name "Travis CI"

# build (CHANGE THIS)
# make

# deploy
cd public
git init
touch index.html
git add index.html
git commit -a -m "Deploy to Github Pages"
#git push --force --quiet "https://${GITHUB_TOKEN}@$github.com/${GITHUB_REPO}.git" master:gh-pages > /dev/null 2>&1
# git push --force --quiet $FULL_REPO master:gh-pages
#git push --force "https://${GITHUB_TOKEN}@$github.com/${GITHUB_REPO}.git" master:gh-pages #> /dev/null 2>&1
git push --force "https://${GITHUB_TOKEN}@$github/${GITHUB_REPO}.git" master:gh-pages #> /dev/null 2>&1

