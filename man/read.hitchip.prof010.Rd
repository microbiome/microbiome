\name{read.hitchip.prof010}
\alias{read.hitchip.prof010}
\title{Description: read hitchip prof 010}
\usage{
  read.hitchip.prof010(path = "../data/", dataID = "")
}
\arguments{
  \item{path}{TBA}

  \item{dataID}{TBA

  Returns:}
}
\value{
  d
}
\description{
  Arguments:
}
\examples{
# usage example:
# d = read.hitchip.prof010("my.file.path")
# d = preproc.Hitchip(d)
# mapply(doPCA,d,names(d))
# where doPCA is a function, for example
# doPCA=function(Data,filename){
#   pdf(paste(filename,"PCAplot.pdf",sep=""))
#   # screeplot:
#   plot(prcomp(t(Data)))
#   dev.off()
#}
#t.test example
# c.i a vector with class information
#t.test.res=apply(d$l2.log10.ave,1,function(x,y) t.test(x[y==1],x[y==2]),y=c.i)
#p.val.vec=sapply(t.test.res,function(x) x$p.value)
#estimate.vec=sapply(t.test.res,function(x) x$estimate)
#qvalue(p.val.vec)$qvalue
#which(p.val.vec<0.05)
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

