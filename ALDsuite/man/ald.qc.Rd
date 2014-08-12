\name{ald.qc}
\alias{ald.qc}

\title{Quality Control Checks}
\description{Checks the quality and the formatting of genetic data.}
\usage{
ald.qc(Pm.prior, haps, geno, gender, chr, pos, burn, iter, every, indiv.id, marker.id, pop.id, lambda, tau, omega, rand.seed, cores, cl, dev, verbose, debug, sex.chr, indiv.geno.check, marker.geno.check, hwe.thresh, bad.indiv, bad.marker, male.het.X, fast, A0, Ak)
}

\arguments{
  \item{Pm.prior}{Prior ancestral information as returned by setup.prior().}
  \item{haps}{Integer matrix with 0/1/2 haplotypes (one row per
    individual, one column per marker).}
  \item{geno}{Integer matrix with 0/1/2 genotypes (one row per
    individual, one column per marker).}
  \item{gender}{Character vector (M/F) - corresponding to the rows of
    geno.}
  \item{chr}{Vector of chromosome number for each marker - corresponding
    to the columns of geno.}
  \item{pos}{Vector of chromosome positions, in cM, for each marker -
    corresponding to chr.}
  \item{burn}{Number of burn-in iterations.}
  \item{iter}{Number of follow-on iterations.}
  \item{every}{Save the state after "every" follow-on iterations.}
  \item{indiv.id}{Character vector of individual id's corresponding to
    the rows of geno.}
  \item{marker.id}{Character vector of marker id's corresponding to the
    columns of geno.}
  \item{pop.id}{Character vector of population id's.}
  \item{lambda}{Estimated sample-wide mean number of generations since
    admixture.}
  \item{tau}{Population dispersion parameter - generally about 300 for
    African Americans.}
  \item{omega}{The parameter describing the dirichlet distributed sample
    prior for admixture (see Details).}
  \item{rand.seed}{Random seed for reproducibility. If NULL, a random
    seed is chosen for you.}
  \item{cores}{Integer number of computing cores to use when
    parallelizing. A value of 0 indicates the program is to be run serially.}
  \item{cl}{A cluster object returned from
    \code{\link{makeCluster}}. When NULL, an appropriate object will be
    generated conditional on the number of cores specified.}
  \item{dev}{Logical value indicating wether the current working version
    or the development version should be run. This argument is generally
    ignored.}
  \item{verbose}{Logical flag used to return information on the
    follow-on iterations and details of data checks.}
  \item{debug}{Logical flag used to return extra output when running in
    verbose mode, including burn-in information.}
  \item{sex.chr}{A scalar value indicating the sex chromosome.}
  \item{indiv.geno.check}{The proportion of genotypes allowed to be
    missing in an individual.}
  \item{marker.geno.check}{The proportion of genotypes allowed to be
    missing for a marker.}
  \item{hwe.thresh}{Hardey Weinberg Equilibrium p-value threshold for
    flagging markers.}
  \item{bad.indiv}{Numeric vector designating rows of geno, elements of
    indiv.id, etc... to drop from the analysis.}
  \item{bad.marker}{Numeric vector designating columns of geno, elements
    of marker.id, etc... to drop from the analysis.}
  \item{male.het.X}{The number of allowable heterozygous loci to allow
    in males without throwing an exception.}
  \item{fast}{Run faster implementation. WARNING: the faster
    implementation is less reproducible.}
  \item{A0}{A matrix with estimated individual admixture proportions
    from each population.}
  \item{Ak}{An array with estimated admixture proportions for each
    individual's parents.}
}

\references{}
\author{Emily Mulhern}


\seealso{see admixture}


\keyword{debugging}

