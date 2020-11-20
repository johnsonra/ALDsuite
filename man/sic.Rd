\name{sic}

\alias{sic}

\title{Calculate Shannon Information Content}

\description{Calculates Shannon Information content.}

\usage{sic(fa, na, fe, ne, admix = 0.2, type = "sic")}

\arguments{
  \item{fa}{Frequency of the African variant.}
  \item{na}{Number of African samples.}
  \item{fe}{Frequency of the European variant.}
  \item{ne}{Number of European samples.}
  \item{admix}{Admixure proportion in African Americans.}
  \item{type}{Type of measure to return (options are \code{sic},
    \code{In}, and \code{Ia}).}
}

\value{A numeric vector.}

\references{
  Rosenberg et. al. Informativeness of Genetic Markers for
  Inference of Ancestry, American Journal of Human Genetics 73:1402--22,
  2003

  Smith et. al. A High-Density Admixture Map for Disease Gene Discovery
  in African Americans, American Journal of Human Genetics 74:1001--13, 2004
}

\author{Randall Johnson}

\note{I have included \code{admix} so the same measures can be
  calculated for other populations.}

\examples{
sic(.2, 134, .8, 207)
}

\keyword{univar}
