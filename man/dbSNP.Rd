\name{dbSNP}

\alias{dbSNP}

\title{Retrieve polymorphism information from dbSNP}

\description{Retrieve polymorphism information from dbSNP}

\usage{dbSNP(snp, parallel = FALSE, processes = 2, type = 'MPI')}

\arguments{
  \item{snp}{A character vector of polymorphism identifiers. Idnetifier must
    begin with `rs' or `hCV', and end in a number (eg rs333 or
    hCV1944294).}
  \item{parallel}{If true, the \code{snow} package is requred, and a
    parallel version will be run.}
  \item{processes}{The number of processes to spawn if \code{parallel}
    is TRUE (see Note).}
  \item{type}{String indicating the type of parallel machine (ie 'MPI'
    or 'PMV').}
}

\note{Most of the time in this function is usually spent waiting for
  NCBI (rather than parcing the xml code), so spawning more processes
  than you have processors may actually increase speed (emphasis on the
  may). On my MacBook Pro, the optimum number of processes to use seems
  to be about 6 (ie 3 times the number of processors). While this is
  dependent on more factors than the number of processors (including how
  busy the NCBI servers currently are), the speed up is significant.}

\value{
  A data frame with the following components:
  \item{organism}{Spiecies}
  \item{rsnumber}{RS number of the polymorphism}
  \item{type}{Type of polymorphism (eg exact (ie SNP), insertion)}
  \item{het}{Heterozygosity estimate}
  \item{het.se}{Standard error of the heterozygosity estimate}
  \item{tag5}{Tagging sequence on the 5' end}
  \item{polymorph}{Polymrphic region (eg -/GTCAGTATCAATTC, A/G)}
  \item{tag3}{Tagging sequence on the 3' end}
  \item{genomeBuild}{Genome build number}
  \item{dbSnpBuild}{dbSNP build number}
  \item{chr}{Chromosome}
  \item{phys.pos}{Physical position of the polymorphism}
}

\author{Randy Johnson}

% \examples{

% x <- dbSNP(c("hCV1944294", 'rs1057079', 'rs333'))

% }

\keyword{database}