\name{show.chain}
\alias{show.chain}

\title{Plot a parameter from the MCMC chain}

\description{This function plots a specific parameter over all iterations of
  the MCMC chain, as recorded by admixture() when the debug option is set.}

\usage{
show.chain(ncpath, param, ...)
}

\arguments{
  \item{ncpath}{A character string indicating the path to a NetCDF file created by admixture()}
  \item{param}{A character string indicating the name of the parameter to be shown}
  \item{\dots}{Any other graphical parameters the user wishes to have passed to plot()}
}

\details{The following parameters are available
  for inspection (see admixture() for more details): A, Ak, AX, lambda,
  lambdaX, P, Pm, omega, omegaX, alpha, alphaX, tau, sigmaA, sigmaL, PCRbetas.
}

\value{The output of this function is a pdf file named after the
  parameter specified by param.
}

\author{
Randall Johnson
}

\examples{
## Not run: show.chain('mald.debug.nc', 'tau')
## Not run: show.chain('mald.debug.nc', 'tau', col = 'blue', pch = 18)
}

\keyword{hplot}