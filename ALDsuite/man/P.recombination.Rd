\name{P.recombination}
\alias{P.recombination}
\alias{P.0crossovers}
\alias{P.1plus.crossover.gamma.gamma.prev}
\alias{P.1plus.crossover.G}
\alias{qodd.crossovers}

\title{Probability of recombination}
\description{
Calculation of the probability of recombination given the number of
generations since admixture and the genetic distance between the markers
of interest, and the probability of no crossovers between two loci.
}

\usage{
P.recombination(x, lambda, distance)
P.0crossovers(lambda, distance)
P.1plus.crossover.gamma.gamma.prev(gamma, gamma.prev, A, lambda, distance)
P.1plus.crossover.G(G, A, lambda, distance, m)
qodd.crossovers(p.crossovers, lambda, distances)
}

\arguments{
  \item{x}{A logical value (TRUE/FALSE) to calculate the probability for
    (e.g. for probability of no recombination, x = FALSE).}
  \item{lambda}{Numeric vector of the numbers of generations since admixture.}
  \item{distance}{Numeric vector of the genetic distances (cM / 100)
    between markers of interest.}
  \item{distances}{Numeric vector of the genetic distances (cM / 100)
    between markers of interest.}
  \item{gamma}{A numeric vector of ancestral population of origin for
    the marker of interest.}
  \item{gamma.prev}{A numeric vector of ancestral population of origin for
    the marker previous to the marker of interest in the MCMC.}
  \item{A}{A numeric vector of global ancestry estimates for one
    individual. A is assumed to sum to 1.}
  \item{G}{An array as returned by sample.G()}
  \item{m}{The number of markers. Equal to dim(G)[1] and
    length(distance).}
  \item{p.crossovers}{Probability to calculate the quantile for.}
}

\details{
  A single value for x is expected. The vector lengths of lambda
  and distance are assumed to either be equal or equal to 1. Lambda or
  distance may be of different lengths provided one or the other is of
  length 1.

  For P.1plus.crossover.gamma.gamma.prev, the length of gamma vectors is
  assumed to be equal to that of lambda and distance. All data are
  assumed come from one individual, who's global ancestry estimates are
  described by A.
}

\value{A vector of probabilities.}

%\references{}

\author{Randall Johnson}

\note{
  A distinction is made here between observing a crossover, which we
  term a recombination and the number of actual crossovers
  between two loci. An even number of crossovers between two loci will
  not result in an observed recombination event, while an odd number of
  crossovers will.

  These functions are not intended for use by the end user, but is
  included here to facilitate ease of development and parallel
  computation.

  P.recombination called by: \code{\link{P.gammas.a}} and
  \code{\link{P.gamma.aX}}.

  P.0crossovers called by: \code{\link{P.1plus.crossover.gamma.gamma.prev}}.

  P.1plus.crossover.gamma.gamm.prev called by: \code{\link{updt.lambda}}.
}

\examples{
P.recombination(TRUE, 6, c(0.1, 0.3, 0.98))
P.recombination(FALSE, 3:6, 0.18)
}

\keyword{distribution}
