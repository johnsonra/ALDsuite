\name{HapMap}

\alias{hapmap36}

\docType{data}

\title{HapMap genetic maps for build 36}

\description{Location map of the human genome with genetic and physical
  positions.}

\usage{data(hapmap36)}

\format{
  hapmap36 is a data frame with 3,281,323 observations on the following 2 variables.
  \describe{
    \item{\code{phys.pos}}{Physical position in the human genome}
    \item{\code{cM}}{Genetic position (averaged over both sexes)}
  }
}

\source{http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/}

\references{
  HapMap Phase III

  HapMap Phase II

  The International HapMap Consortium. The International HapMap Project. Nature 426, 789-796 (2003)}

\examples{
data(hapmap36)
}

\keyword{datasets}
