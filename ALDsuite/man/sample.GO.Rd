\name{sample.G}
\alias{sample.G}
\alias{sample.G.par}
\alias{sample.O}
\alias{updt.G0}

\title{Sample ancestral states.}

\description{
  Functions and helper functions for sampling ancestral states
  (sample.G) and the ancestral state associated with each allele
  (sample.O).
}

\usage{
sample.G(gammas.joint, pops, chr, gender, sex.chr, cl)
sample.G.par(j, gammas.joint, chr, gender, sex.chr)

sample.O(geno, P, G)

updt.G0(G)
}

\arguments{
  \item{gammas.joint}{An array of gamma values returned from
    \code{\link{sample.joint.gammas}}.}
  \item{pops}{A vector of populations under study
    (i.e. \code{dimnames(P)[[2]]}).}
  \item{chr}{A vector of chromosome numbers for each locus.}
  \item{gender}{A vector of gender specifications ("M" or "F") for each
    individual.}
  \item{sex.chr}{A scalar value indicating which chromosome is the sex
    chromosome (e.g. 23, "X", ...).}
  \item{cl}{An object of class 'c("SOCKcluster", "cluster")' as returned
    by \code{\link{makeCluster}}. Set to NULL to run serially.}
  \item{j}{A scalar index indicating which marker to sample.}
  \item{geno}{A matrix of genotypes as defined for
    \code{\link{sample.joint.gammas}}.}
  \item{P}{A matrix of estimated modern allele frequencies from each
    ancestral population as defined for \code{\link{sample.joint.gammas}}.}
  \item{G}{An array of sampled ancestral states as returned by
    \code{\link{sample.G}}.}
}

\value{
  An array of values...fill this in!
}

%\references{}

\author{Randall Johnson}

\note{
  These are not intended for use by the end user, but is included here
  to facilitate ease of development and parallel computation.
}

%\examples{}

\keyword{distribution}
