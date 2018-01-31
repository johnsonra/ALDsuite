\name{plot.amap}

\alias{plot.amap}

\title{Plot amap}

\description{Plot an amap object}

\usage{
\method{plot}{amap}(x, cM = FALSE, add.legend = TRUE, chromosome, ci =
0.95, ci.col = "grey75", plot.type = "genome", change.par = TRUE, ...)
}

\arguments{
  \item{x}{An object of class amap}
  \item{cM}{logical indicating x axis scale. Default is to print on Mb
    scale (cM is not currently implemented)}
  \item{add.legend}{logical indicating wether a legend should be added
    to the plot}
  \item{chromosome}{-}
  \item{ci}{-}
  \item{ci.col}{-}
  \item{plot.type}{-}
  \item{change.par}{Logical value indicating wether arguments to par()
    should be changed}
  \item{...}{Other graphical parameters (none implemented)}
}

\author{Randall Johnson <rjohnson@ncifcrf.gov>}

\examples{
  data(anc.ex1)
  plot(anc.ex1)
}

\keyword{ hplot }% at least one, from doc/KEYWORDS
