\name{rdirichlet}
\alias{rdirichlet}

\title{Random dirichlet sample}

\description{Returns a dirichlet sample in matrix form.}

\usage{
rdirichlet(n, alpha)
}

\arguments{
  \item{n}{The number of samples.}
  \item{alpha}{Shape parameters for the underlying gamma distribution.}
}

\details{This is a simple wrapper function that samples from the
  appropriate gamma distributions and normalizes the results.}

\value{A matrix is returned with n rows and length(alpha) columns. Each
  row corresponds to an independent sample.}

\author{Randall Johnson}

\examples{
rdirichlet(n = 3, alpha = c(3, 5, 1))
}

\keyword{distributions}
