\name{log likelihoods}

\alias{logP.lambdas.hyper}
\alias{logP.P.tau.Pm}
\alias{logP.Ak.omega}
%\alias{logP.omega.Ak}
\alias{logP.AX.omega}
%\alias{logP.omega.AX}
\alias{logP.gammas.geno.P}
\alias{logP.gammas.geno.P.j}
\alias{logP.lambda.gammas}
\alias{logP.delta.sigma}

\title{Log likelihood functions}

\description{
  Log likelihood functions for use in updates.
}

\usage{
logP.lambdas.hyper(lambdas, alpha, beta, epsilon = 0.01)

logP.P.tau.Pm(P, tau, Pm)

logP.Ak.omega(Aks, omega, ratio = TRUE)
%logP.omega.Ak(Aks, omega)

logP.AX.omega(AX, Aks, omegaX, ratio = TRUE)
%logP.omega.AX(AX, Aks, omegaX)

logP.gammas.geno.P(gammas, geno, P)
logP.gammas.geno.P.j(j, gammas, geno, P)

logP.lambda.gammas(gammas, lambda, distances)

logP.delta.sigma(deltas, sigma, lims, dev, epsilon = 0.01)
}

\arguments{
  \item{lambdas}{Vector of lambda estimates for each individual.}
  \item{lambda}{Numeric value for lambda.}
  \item{alpha}{Scale parameter for prior on lambda.}
  \item{beta}{Rate parameter for prior on lambda.}
  \item{epsilon}{A minimum value for lambda/delta to avoid returning
    infinite log likelihoods. Should be something close to 0.}
  \item{P}{Vector or matrix of current estimates of modern allele
    frequencies from each ancestral population.}
  \item{tau}{Hyper parameter for prior on P.}
  \item{Pm}{...}
  \item{Aks}{Matrix of local ancestry estimates.}
  \item{AX}{Matrix of local ancestry estimates, specific to the sex
    chromosome.}
  \item{omega}{Hyper parameter for prior on ancestry estimates.}
  \item{omegaX}{Hyper parameter for prior on ancestry estimates,
    specific to the sex chromosome.}
  \item{gammas}{Array of ancestral state probabilities.}
  \item{geno}{Matrix of 0/1/2 genotypes.}
  \item{j}{Integer referring to the marker of interest.}
  \item{distances}{Vector of genetic distances between markers.}
  \item{ratio}{Allows for cancelation of terms when calculating a log
    likelihood ratio, resulting in a quicker evaluation time.}
  \item{deltas}{<fill in>}
  \item{sigma}{<fill in>}
  \item{lims}{Limits on the distribution of deltas. The function returns
    -Inf when delta is outside the limits ... i.e. log(0)}
  \item{dev}{Logical value indicating development code should be
    run. This is generally ignored.}
}

\details{
  Discuss difference in Ps and double check on Pm.

  Figure out dimensions expected for each of these (should be fairly
  consistent) and how they are indexed.
}

\value{
  A log likelihood statistic, or a vector of log likelihood statistics.
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
