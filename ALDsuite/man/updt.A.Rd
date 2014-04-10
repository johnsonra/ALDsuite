\name{updt.A}
\alias{updt.A}
\alias{updt.omega}
\alias{updt.AX}
\alias{updt.omegaX}
\alias{updt.sigma}

\title{Update Ancestry Estimates}

\description{
  Global ancestry estimates are updated using updt.A, ancestral state
  probabilities for an individual are calculated using updt.Aj.
}

\usage{
updt.A(A0, Ak, Ajs, omega, sigma, epsilon, dev)
updt.AX(A0, AX, Ak, Ajs, omegaX, sigma, W, dev)

updt.omega(omega, A0, burn = 50, iter = 100, every = 10, debug = FALSE)
updt.omegaX(omegaX, Ak, AX, M, burn = 50, iter = 100, every = 10, debug = FALSE)

updt.sigma(sigma, deltas, lims, dev, burn = 50, iter = 100, isA = TRUE)
}

\arguments{
  \item{A0}{Array of global ancestry estimates for each individual (see
    Details).}
  \item{Ak}{Array of global ancestry state estimates for each
    parent (see Details).}
  \item{AX}{Array of global ancestry state estimates for each
    parent for the sex chromosome (see Details).}
  \item{Ajs}{Array of ancestral state probabilities for each individual
    (see Details).}
  \item{epsilon}{A minimum value for A before A is considered equal to 1
    or 0 (i.e. no admixture for the parent or individual). Should be
    something close to 0.}
  \item{omega}{Vector of values defining the Dirichlet prior for each
    individual's global ancestry (see Details).}
  \item{omegaX}{Vector of values defining the Dirichlet prior for each
    individual's global ancestry on the sex chromosome (see Details).}
  \item{W}{Logical vector of the same length as the number of
    individuals in Ak, indicating if the individual is female.}
  \item{M}{Logical vector of the same length as the number of
    individuals in Ak, indicating if the individual is male.}
  \item{burn}{Number of burn-in iterations for update.}
  \item{iter}{Number of follow-on iterations for update.}
  \item{every}{Sample the chain after 'every' cycles.}
  \item{debug}{Logical indicating how much information to return.}
  \item{sigma}{Hyper parameter indicating the spread of parents'
    estimates around the indivdiual's global ancestry (i.e. how different
    are the parents' global ancestries from each other in the population
    under study).}
  \item{deltas}{Inferred difference in global ancestry of the parents of
    an individual.}
  \item{isA}{Logical: TRUE indicates sigma is being updated for A, FALSE
    indicates it is being updated for lambda.}
  \item{lims}{The limits on the possible range of deltas (e.g. delta
    could not be 0.5 when the individual's estimate is 0.1).}
  \item{dev}{Logical value indicating development code should be
    run. This is generally ignored.}
}

\details{
  A0 has 2 dimensions: 1 - Individual, 2 - Population.

  Ak has 3 dimensions: 1 - Individual, 2 - Parent, 3 - Population.

  Ajs has 3 dimensions: 1 - Individual, 2 - Marker, 3 - Population.

  distance is a vector of genetic distances between markers, where the
  first marker of each chromosome has a missing value.

  omega has length equal to the number of populations under study
  (i.e. dim(A)[2] == dim(Ajs)[3] == length(omega)).
}

\value{
  An array of updated values is returned with the same dimensions as A
  and omega, respectively.
}

%\references{}

\author{Randall Johnson}

\note{
  These functions are not intended for use by the end user, but are
  included here to facilitate ease of development and parallel
  computation.
}

%\examples{}

\keyword{distribution}