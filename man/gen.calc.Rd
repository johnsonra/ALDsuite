\name{gen.calc}

\alias{gen.calc}
\alias{phys.calc}

\title{Calculate Physical/Genetic Map Position}

\description{Calculate physical/geneetic map position based on the Matise
  et. al. genetic map (build 36).}

\usage{gen.calc(chrom, pos, n.extrap = 4)
       phys.calc(chrom, pos, n.extrap = 4)}

\arguments{
  \item{chrom}{vector of chromosome numbers}
  \item{pos}{vector of (build 36) physical/genetic positions
  (gen.calc/phys.calc respectively)}
  \item{n.extrap}{The number of neighboring points to use when extrapolating}
}

\details{This function expects \code{chrom} to be numeric (ie the X
  chromosome needs to be represented by 23), and uses the \code{rutgers}
  data set included in this package to interpolate the physical/genetic map
  position based on the genetic/physical map position.}

\value{
  A data frame with 3 variables:
  \item{chr}{Chromosome number}
  \item{phys.pos}{Physical map position}
  \item{gen.pos}{Genetic map position}
}

\references{See the following paper on the development of the rutgers
  map: Matise et.al, {American Journal of Human Genetics},
  2003, 73(2):271-284}

\author{Randy Johnson}

\seealso{\code{\link{rutgers}}}

\keyword{database}
