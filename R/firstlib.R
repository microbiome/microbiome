# Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi Contact:
# <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or modify it
# under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

#' @import knitr
.onAttach <- function(lib, pkg) {
    packageStartupMessage("\nmicrobiome R package (microbiome.github.com)\n\n
        Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi \n
        <microbiome-admin@googlegroups.com>\n")
    
} 
