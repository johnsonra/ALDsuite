\name{read.amap}

\alias{read.amap}

\title{Read ANCESTRYMAP Output}

\description{Read ANCESTRYMAP output files into R}

\usage{read.amap(file.name, remove.files=TRUE)}

\arguments{
  \item{file.name}{Character: Name of output file}
  \item{remove.files}{Logical: Remove temporary files when finished}
}

\value{
  \item{An object of class amap}{List of 7}
}

\references{See Patterson et al, Methods for High-Density Admixture Mapping of Disease Genes, American Journal of Human Genetics
  (74:5, pp 979-1000)}

\author{Randy Johnson}

\seealso{See also \code{\link{amap}} for setting up and running
  ANCESTRYMAP. Web pages with more information are also available at}

\keyword{manip}
\keyword{IO}
