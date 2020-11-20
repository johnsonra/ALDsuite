\name{updt.lambda}
\alias{updt.lambda}
\alias{updt.alpha.beta}
\alias{sample.lambda}
\alias{sample.lambda.k}

\title{Updates for lambda}

\description{
  Update estimates for the number of generations since admixture and
  other associated functions.
}

\usage{
updt.lambda(G, G0, A0, Ak, lambda, sigma, distances, alpha, beta, chr, P.Aeq1, dev, W = NULL)

updt.alpha.beta(alpha, beta, lambdas, step = 1, burn = 50, iter = 100,
every = 10, debug = FALSE)

sample.lambda(G, A, distances, chr, alpha, beta, P.Aeq1, dev)
sample.lambda.k(G0, G, A, distances, chr, alpha, beta, P.Aeq1, dev)
}

\arguments{
  \item{G}{Sampled ancestral states as returned by sample.G().}
  \item{G0}{Sampled ancestral states as returned by sample.G0().}
  \item{A0}{Array of global ancestry estimates for each individual (see
    Details).}
  \item{lambdas}{Matrix of lambda values (see Details).}
  \item{lambda}{Matrix of lambda values (see Details).}
  \item{sigma}{Hyper parameter indicating the spread of parents'
    estimates around the indivdiual's value for lambda (i.e. how different
    are the parents' lambdas from each other in the population
    under study).}
  \item{distances}{Vector of distances between each marker.}
  \item{alpha}{Shape parameter for prior for lambda.}
  \item{beta}{Rate parameter for prior for lambda.}
  \item{Ak}{Array of global ancestry state estimates for each
    parent (see Details).}
  \item{chr}{A vector of chromosome numbers for each marker.}
  \item{A}{Vector gloabal ancestry estimates for one parent (see
    Details).}
  \item{P.Aeq1}{Array of probabilities that a parent has no admixture
    (see Details).}
  \item{step}{Parameter for scaling the step size between guesses.}
  \item{burn}{Number of burn-in iterations for update.}
  \item{iter}{Number of follow-on iterations for update.}
  \item{every}{Sample the chain after 'every' cycles.}
  \item{debug}{Logical indicating how much information to return.}
  \item{W}{Logical vector of the same length as the number of
    individuals in G, indicating the individual is female (FALSE
    indicates the individual is male).}
  \item{dev}{Logical value indicating development code should be
    run. This is generally ignored.}
}

\details{
  A0 has 2 dimensions: 1 - Individual, 2 - Population.

  lambdas is a matrix with the following dimensions: 1 - Individual, 2 -
  Parent.

  Ak is an array with the following dimensions: 1 - Individual, 2 -
  Parent, 3 - Population.

  A is a vector of length equal to the number of populations under
  study.

  P.Aeq1 is an array of probabilities that a parent has no admixture
  from a given population. It is returned from updt.A and has the
  following dimensions: 1 - Individual, 2 - Parent, 3 - Population. Note
  that the value passed to sample.lambda.k will be missing the second
  dimension (i.e. it contains values for only one parent).
}

\value{
  update.lambda.par returns an array of updated values with the same
  dimensions as lambdas.

  update.lambda returns an updated value for lambda.

  update.alpha.beta returns a named list with a updated values for alpha
  and beta.
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
