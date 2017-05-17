#!/bin/bash

SOURCE_BRANCH="devel"
TARGET_BRANCH="master"

# Save some useful information
#SHA=`git rev-parse --verify HEAD`
GH_REPO="@github.com/microbiome/microbiome.git"
FULL_REPO="https://$GH_TOKEN$GH_REPO"

# config
git config --global user.email "nobody@nobody.org"
git config --global user.name "Travis CI"

# Add to git and deploy
git fetch origin master
git merge master
git commit -a -m "Merge devel with master"
git push --quiet $FULL_REPO $SOURCE_BRANCH:$TARGET_BRANCH

