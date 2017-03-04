#!/bin/bash

set -o errexit -o nounset

if [ "$TRAVIS_BRANCH" != "master" ]
then
  echo "This commit was made against the $TRAVIS_BRANCH and not the master! No deploy!"
  exit 0
fi

rev=$(git rev-parse --short HEAD)

cd stage/_book

git init
git config user.name "antagomir"
git config user.email "leo.lahti@iki.fi"

git remote add upstream "https://$GH_TOKEN@github.com/microbiome/microbiome.git"
git fetch upstream
git reset upstream/gh-pages

echo "microbiome.com" > index.html

touch .

git add -A .
git commit -m "rebuild pages at ${rev}"
git push -q upstream HEAD:gh-pages

