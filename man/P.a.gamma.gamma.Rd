\name{P.a.gamma.gamma}
\alias{P.a.gamma.gamma}

\title{Probability of a genotype}

\description{This function calculates the probability of a
  genotype, conditional on both parent's ancestral origin at the locus.}

\usage{
P.a.gamma.gamma(x, gamma, gamma.alt, P)
}

\arguments{
  \item{x}{A vector of 0/1/2 genotypes.}
  \item{gamma}{A vector of ancestral origins for one parent.}
  \item{gamma.alt}{A vector of ancestral origins for the other parent.}
  \item{P}{A vector of modern allele frequencies in loci from each
    ancestral population.}
}

\details{
  The vector lengths for x, gamma and gamma.alt should all be equal, and
  values for gammas should be in 1:length(P).

  Genotypes in x should all be from the same marker.
}

\value{A vector of probabilities.}

%\references{}

\author{Randall Johnson}

\note{
  These are not intended for use by the end user, but is included here
  to facilitate ease of development and parallel computation.

  P.a.gamma.gamma called by: \code{\link{P.gammas.a}}.
}

%\examples{}

\keyword{distribution}
