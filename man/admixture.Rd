\name{admixture}
\alias{admixture}

\title{Inferrence of individual ancestry}

\description{
  Individual ancestry is inferred at each marker along the genome.
}

\usage{
admixture(geno, gender, chr, pos = NULL, build = 37, d = NULL, Pm.prior,
burn = 100, iter = 200, every = 2, indiv.id = NULL,
marker.id = NULL, pop.id = NULL, lambda = 9, tau = 300, omega = NULL,
rand.seed = NULL, cores = detectCores(), cl = NULL, dev = FALSE,
verbose = TRUE, debug = FALSE, sex.chr = 23, indiv.geno.check = 0.98,
marker.geno.check = 0.98, hwe.thresh = 1e-4, bad.indiv = NULL,
bad.marker = NULL, male.het.X = 1, fast = FALSE, coolness = 1)
}

\arguments{
  \item{geno}{Integer matrix with 0/1/2 genotypes (one row per
    individual, one column per marker).}
  \item{gender}{Character vector (M/F) - corresponding to the rows of
    geno.}
  \item{chr}{Vector of chromosome number for each marker - corresponding
    to the columns of geno.}
  \item{pos}{Vector of chromosome positions, in cM, for each marker -
    corresponding to chr.}
  \item{build}{Build of human genome used for pos (supported builds are
    36 and 37)}
  \item{d}{Numeric vector of genetic distances (in cM) between each
    marker and the previous marker (first marker on a chromosome is
    NA). Either pos or d must be specified.}
  \item{Pm.prior}{Prior ancestral information as returned by setup.prior().}
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
  \item{coolness}{A numeric value in [0,1] restricting how quickly P can
    change. This is generally not needed.}
}

\details{
  omega: This can either be a single, positive numeric value (indicating
  a flat prior) or a vector of the same length as pop.id.
}

\value{
  A list is returned with a varying number of elements, depending upon
  values for the verbose and debug flags. The first element, "final", is
  alway included.

  final:
  \item{A0}{A matrix with estimated individual admixture proportions
    from each population.}
  \item{Ak}{An array with estimated admixture proportions for each
    individual's parents.}
  \item{Aj}{An array of locus specific admixture state probabilities for
    each individual.}
  \item{P}{A matrix of ancestry-specific allele frequency estimates in
    the current population.}
  \item{omega}{A vector of parameters describing the prior distribution
    of global admixture porportions from each population.}
  \item{lambda}{A matrix of estimated mean number of generations since
    admixture for each individual's parents.}
  \item{alpha}{The first hyper-parameter describing the prior
    distribution of lambda in the sampled population.}
  \item{beta}{The second hyper-parameter describing the prior
    distribution of lambda in the sampled population.}
  \item{tau}{A vector of parameters describing the prior distribution of
    ancestry-specific allele frequencies in the current population (one
    value of tau for each ancestral population).}
  \item{omegaX}{A vector of parameters describing the prior distribution
    of X chromosome admxture proportions from each population.}
  \item{lambdaX}{A matrix of estimated effective mean number of
    generations since admixture for the X chromosome.}
  \item{alphaX}{The first hyper-parameter describing the prior
    distribution of lambdaX in the sampled population.}
  \item{betaX}{The first hyper-parameter describing the prior
    distribution of lambdaX in the sampled population.}

  Other items possibly included in the list are:
  \item{burn}{A list with items similar to final, but with samples taken
    during the burnin phase.}
  \item{iter}{A list with items similar to final, but with samples taken
    after the burnin phase.}
  \item{sex.check}{A named vector with statistics for individuals
    flagged during the sex check.}
  \item{indiv.check}{A named vector with statistics for individuals
    flagged during the analysis of missing data.}
  \item{marker.check}{A named vector with statistics for makers flagged
    during the analysis of missing data.}
  \item{hwe.check}{A named vector with statistics for markers flagged
    during the HWE check.}
}

\note{
  It is a good idea to include dimnames for your data objects and values
  for indiv.id, marker.id and pop.id. Inclusion of row names, column
  names, etc... will allow the software to double check that everything
  is sorted the way you think it is.

  Due to the way a cluster generates random numbers, it is unwise to
  assume a random seed will produce exactly the same results on clusters
  of varying sizes.
}

%\references{}

\author{Randall Johnson}

\examples{
# See the viniette for this package for examples.
}

\keyword{optimize}
