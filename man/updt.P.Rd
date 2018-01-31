\name{updt.P}
\alias{updt.P}
\alias{updt.tau.Pm}

\title{Updates of Estimated Allele Frequencies}

\description{
  These functions are for updates of estimated allele frequencies in
  ancestral populations, the admixed population, and associated
  parameters.
}

\usage{
updt.P(P, tau, Pm, Os, coolness, dev)
updt.tau.Pm(tau, Pm, Pm.counts, P, burn = 50, iter = 100, every = 10, debug = FALSE)
}

\arguments{
  \item{P}{Matrix of reference allele frequencies (see Details).}
  \item{tau}{Vector of population dispersion factors (one per population).}
  \item{Pm}{Matrix of ancestral reference allele frequencies (see
    Details).}
  \item{Pm.counts}{Allele frequency counts in surrogate ancestral
    populations (see Details).}
  \item{Os}{...sampled alleles with population designation...}
  \item{coolness}{A numeric value in [0,1] restricting how quickly P can
    change (see Details).}
  \item{burn}{Number of burn-in iterations for update.}
  \item{iter}{Number of follow-on iterations for update.}
  \item{every}{Sample the chain after 'every' cycles.}
  \item{debug}{Logical indicating how much information to return.}
  \item{dev}{Logical indicating development code should be run. This is
    generally ignored.}
}

\details{
  P and Pm are matrices with the following dimensions: 1 - Marker, 2 -
  Population.

  Pm.counts is an array with the following dimensions: 1 - Marker, 2 -
  Reference / Variant, 3 - Population.

  The coolness variable was added as an option to slow the rate of
  change of P from one iteration to the next. 0 results in P not being
  updated at all, and 1 results in updating to the new value without
  taking the old value into consideration. Intermediate values result in
  a weigted sum of the old and new values.
}

\value{
  updt.P returns an updated matrix of allele frequencies with the same
  dimensions as P.

  updt.tau.Pm returns a named list:
  \item{Pm}{an updated matrix of allele frequencies with the same
    dimensions as Pm.}
  \item{tau}{a vector of population dispersion factors with the same
    length as tau.}
}

%\references{}

\author{Randall Johnson}

\note{
  These are not intended for use by the end user, but is included here
  to facilitate ease of development and parallel computation.
}

%\examples{}

\keyword{distribution}
