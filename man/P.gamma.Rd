\name{P.gamma}
\alias{P.gamma}
\alias{P.gammas.a}
\alias{P.gammas.marg}
\alias{P.gamma.aX}

\title{Prbability of a set of gammas.}

\description{
  Probability that a locus originated from ancestral population
  x. P.gammas.a is conditional on genotype and calculates the joint
  probabilities, including both parents. P.gammas.marg calculates
  the marginal probabilities for each individual parent, independent of
  the other, based on the joint probabilities calculated using
  P.gammas.a. P.gamma.aX calculated the probabiltiy of gamma on the male
  X chromosome.
}

\usage{
P.gamma(x, lambda, distance, Ak, gamma.prev = rep(NA, dim(Ak)[2]), p.recomb.T, p.recomb.F)

P.gammas.a(combos, a, P, lambda, d, Ak, gammas.prev)

P.gammas.marg(gammas.joint, pops)

P.gamma.aX(x, a, P, lambda, d, Ak, gamma.prev)
}

\arguments{
  \item{x}{An ancestral population code (e.g. 1, 2, ...).}
  \item{lambda}{The numbers of generations since admixture (see Details).}
  \item{distance}{Numeric vector of the genetic distances (cM / 100)
    between markers of interest.}
  \item{Ak}{A matrix of global ancestry estimates for one parent (see
    Details).}
  \item{gammas.prev}{An array of gammas at the previous marker (see Details).}
  \item{gamma.prev}{A matrix of the values for gamma at the previous
    locus for one parent in the MCMC (see Details).}
  \item{p.recomb.T}{Vector of probabilities of recombination.}
  \item{p.recomb.F}{Vector of probabilities of no recombination.}
  \item{a}{A vector of 0/1/2 genotypes.}
  \item{P}{A vector of modern allele frequencies in loci from each
    ancestral population.}
  \item{gammas.joint}{Output from \code{\link{sample.joint.gammas}}.}
  \item{pops}{Populations under study (i.e. \code{dimnames(P)[[2]]} used
    in \code{\link{sample.joint.gammas}}).}
  \item{combos}{Possible ancestral state combinations inherited from
    parents (e.g. \code{c('g11', 'g12', 'g21', 'g22')}).}
  \item{d}{Scalar value of the genetic distance to the previous marker.}
}

\details{
  A single value for x is expected and must fall in 1:dim(Ak)[2]. The
  vector lengths of lambda, distance, p.recomb.T and p.recomb.F as well
  as the number of rows of Ak and gamma.prev should all be equal. The
  numbers of columns of Ak and gamma.prev should be equal.

  dim(Ak)[1] and dim(gamma.prev)[1] are indexed by individual, while
  dim(Ak)[2] and dim(gamma.prev)[2] are indexed by population.

  Ak and gamma.prev may have only one row. In this case vectors of the
  appropriate length may be substituted for matricies, and all other
  values should be scalar values.

  lambda is a vector for a single parent in the case of P.gamma.a
  and P.gamma, but is a two column matrix for P.gammas.a with each
  column representing a parent. In other words, the variable for
  P.gammas.a could be defined with \code{lambda <- cbind(lambda.parent1,
    lambda.parent2)}.

  gammas.prev is an array containing ancestry probabilities or samples
  with the following dimensions: Dim 1 - Individual, Dim 2 - Combos.
}

\value{
  P.gamma and P.gamma.aX each return a vector of probabilities.

  An array of ancestral state likelihoods (for each individual's
  parents, for the marker of interest, in each population under study)
  is returned from P.gammas.a with the following dimensions:
  \item{Dim1}{Individual}
  \item{Dim2}{Parent}
  \item{Dim3}{Combination of potential parent's ancestral states
    (e.g. "g11", "g12", "g21", "g22").}

  An array of ancestral state likelihoods (for each individual's
  parents, for each marker, in each population under study) is returned
  from P.gammas.marg with the following dimensions:
  \item{Dim 1}{Individual}
  \item{Dim 2}{Marker}
  \item{Dim 3}{Parent}
  \item{Dim 4}{Population}
}

%\references{}

\author{Randall Johnson}

\note{
  These functions are not intended for use by the end user, but is
  included here to facilitate ease of development and parallel
  computation.

  P.gamma called by: P.gamma.aX.

  P.gammas.a called by: \code{\link{sample.joint.gammas.par}}.

  P.gammas.marg called by: \code{\link{admixture}} and
  \code{\link{sample.joint.gammas.par}}.
}

%\examples{}

\keyword{distribution}
