# Development instructions for command line:

## Git basics

1. Clone the [git repository](https://github.com/microbiome/microbiome/):  
<pre><code>git clone git@github.com:microbiome/microbiome.git</pre></code>

1. Pull changes by others:  
<pre><code>git pull</pre></code>

1. Make your changes  

1. Document, build, check and install the modified package to check its validity (see below)

1. Make [pull request](https://github.com/microbiome/microbiome) or push directly to the repository (for pushing you need permissions; please contact the [admins](contact)):  
<pre><code>git add new.files
   git commit -a -m"my changes"
   git push
</pre></code>

See also [source code at GitHub](https://github.com/microbiome).

## Building the R package

Instructions for building the R package tarball, assuming that the package structure is under the "pkg" directory in the working directory.

1. Load the devtools R package (requires >=R-2.15)
<pre><code>library("devtools")</pre></code>

1. Build Roxygen documentation
<pre><code>document("pkg")</pre></code>

1. Validate the package
<pre><code>check("pkg")</pre></code>

1. Build the package
<pre><code>build("pkg")</pre></code>

1. Install the package
<pre><code>install("pkg")</pre></code>

When ready: add, commit and push to git master repository.

