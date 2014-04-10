\name{sample.gammas}
\alias{sample.joint.gammas}
\alias{sample.joint.gammas.par}
\alias{updt.Aj}

\title{Calculation of ancestral state likelihoods}

\description{
  These functions calculate the ancestral state likelihoods as a function
  of other parameters and observed data. sample.joint.gammas does the
  calculation for each parent and updt.Aj gives combined probabilities
  for the individual (parents are combined).
}

\usage{
sample.joint.gammas(geno, P, lambda, lambdaX, d, chr, Ak, AX, cl,
                    gender, sex.chr, dev)
sample.joint.gammas.par(chr.curr, geno, P, lambda, lambdaX, d, chr, Ak,
                        AX, combos, gender, sex.chr, dev)

updt.Aj(gammas.joint)
}

\arguments{
  \item{geno}{Matrix of genotypes (see Details).}
  \item{P}{A matrix of estimated modern allele frequencies from each
    ancestral population (see Details).}
  \item{lambda}{A matrix of lambda values indicating the number of
    generations since admixture for each parent (see Details)}
  \item{lambdaX}{A matrix identical in structure to lambda, but with
    values for the X chromosome (see Details).}
  \item{d}{A vector of genetic distances between markers.}
  \item{chr}{A vector of chromosome numbers. Should be the same length
    as d.}
  \item{Ak}{An array of global ancestry estimates (see Details).}
  \item{AX}{An array identical in structure to Ak, but with values for
    the X chromosome (see Details).}
  \item{cl}{A cluster object as returned by \code{\link{makeCluster}},
    or NULL if it is to be run serially}
  \item{gender}{A character or factor vector indicating gender. Values
    must be "M" or "F".}
  \item{sex.chr}{A numeric for character scalar value indicating which
    is the sex chromosome. If sex chromosome data are not given, any
    value not in chr may be passed.}
%  \item{gammas}{Output from \code{\link{P.gammas.marg}} or NULL.}
  \item{chr.curr}{A scalar value indicating the chromosome to sample.}
  \item{combos}{A vector of all possible ancestral combinations,
    coninciding with the third dimension of the output of
    sample.joint.gammas.}
  \item{gammas.joint}{Output from sample.joint.gammas.}
  \item{dev}{Logical value indicating development code should be
    evaluated.}
}

\details{
  Ak - dim1 = individual, dim2 = parent, dim3 = population

  AX - This is identical in structure to Ak, but the second
  column will be NA for males, since they have only one X chromosome.

  geno -

  P -

  lambda -

  lambdaX - This is identical in structure to lambda, but the second
  column will be NA for males, since they have only one X chromosome.

  d -
}

\value{
  A similar array of ancestral state likelihoods is returned from
  \code{\link{sample.joint.gammas}}, but the 3rd and 4th dimensions are
  combined. In this object, combinations of populations inherited are
  returned instead, similar to the output from \code{\link{updt.Aj}},
  but where the order matters. An example of the dimnames is
  \code{c('g11', 'g12', 'g21', 'g22')}. It has the following dimensions:
  \item{Dim 1}{Individual}
  \item{Dim 2}{Parent}
  \item{Dim 3}{Combination of populations inherited (one from each
    parent)}

  An array of ancestral state likelihoods is returned from
  \code{\link{updt.Aj}}. This differs from the output of
  \code{\link{sample.joint.gammas}} in that the order of populations
  does not matter between parents. An example of the dimnames of the 3rd
  dimension is \code{c('g11', 'g12', 'g22')}. It has the following
  dimensions:
  \item{Dim 1}{Individual}
  \item{Dim 2}{Marker}
  \item{Dim 3}{Combination of populations inherited (one from each
    parent)}
}

%\references{}

\author{Randall Johnson}

\note{
  This is not intended for use by the end user, but is included here
  to facilitate ease of development and parallel computation.
}

%\examples{}

\keyword{distribution}
