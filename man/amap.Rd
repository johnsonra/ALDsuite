\name{amap}

\alias{amap}

\title{Run/Setup ANCESTRYMAP files}

\description{amap() creates parameter files and will run ANCESTRYMAP on
  UNIX / LINUX alikes.}

\usage{
amap(param.file.name = "param", preset = NULL, risk = 2,
       indivname = "indiv.dat", snpname = "snps.dat",
       genotypename = "geno.dat", badsnpname = NA, tlreest = TRUE,
       seed = NA, splittau = TRUE, fancyxtheta = TRUE, output = NA,
       trashdir = NA, checkit = FALSE, details = FALSE, numburn = 1,
       numiters = 5, emiter = 10, dotoysim = FALSE, cleanit = TRUE,
       reestiter = 5, thetafilename = NA, lambdafilename = NA,
       freqfilename = NA, indoutfilename = NA, snpoutfilename = NA,
       ethnicfilename = "ethnic.out", noxdata = FALSE, fakespacing = NA,
       thxpars = NA, thpars = NA, lampars = NA, lamxpars = NA,
       markersim = NA, simnumindivs = NA, risksim = NA, tauscal = NA,
       wrisk = NA, lrisk = NA, controlrisk = NA, risk2 = NA, taulsdev = NA,
       taulmean = NA, allmale = NA, allcases = NA, usecontrols = NA,
       pbfmodern = NA, pbxname = NA, genotoyoutfilename = NA,
       indtoyoutfilename = NA, unknowngender = NA, ethnic.control = NA,
       parameter.name = NULL, value = NULL, run = TRUE, in.prefix = NULL,
       out.prefix = NULL, remove.files = TRUE)
}

\arguments{
  \item{param.file.name}{Character string: Parameter file name.}
  \item{preset}{Numeric: Defaults for 3 analyses described in the ANCESTRYMAP
    software documentation. Valid values are: -1 - simulations, 0 -
    input file checks, 1 - data checks, 2 - full run, NULL - All user defined. See
    Details.}
  \item{risk}{Numeric: Vector of risks to model (of length 1 or
    more).}
% ANCESTRYMAP Input File Names
  \item{indivname}{Character: Individual data file name.}
  \item{snpname}{Character: Marker information file name.}
  \item{genotypename}{Character: Individual genotype data file name.}
  \item{badsnpname}{Characer: File name with a list of bad
    SNPs. Default: missing}
% ANCESTRYMAP Main Parameters
  \item{tlreest}{Logical: Always set to TRUE, see ANCESTRYMAP software
    documentation.}
  \item{seed}{Integer: Seed for random number generator; if NA, a random
    seed is chosen.}
  \item{splittau}{Logical: ?}
  \item{fancyxtheta}{Logical: ?}
  \item{output}{Character: Output file for parameter values at every
    iteration. Default: missing}
  \item{trashdir}{Character: Use only in \code{chekit} mode: directory
    to store HMM output. Default: missing}
  \item{checkit}{Logical: Run extra checks.}
  \item{details}{Logical: Generate aditional output.}
  \item{numburn}{Integer: Number of burn-in iterations.}
  \item{numiters}{Integer: Number of iterations.}
  \item{emiter}{Integer: EM iterations.}
  \item{dotoysim}{Logical: Run simulations when TRUE.}
  \item{cleanit}{Logical: ?}
  \item{reestiter}{Integer: Controls the number of iterations inside ANCESTRYMAP for allele frequency sampling.}
% ANCESTRYMAP output file names
  \item{thetafilename}{Character: Ancestry output file for all
  individuals. Default: "theta.out" in non-\code{checkit} mode; missing otherwise.}
  \item{lambdafilename}{Character: Lambda output file for all
  individuals. Default: "lambda.out" in non-code{checkit} mode; missing otherwise.}
  \item{freqfilename}{Character: Allele frequency output file. Default:
  "freq.out" in non-code{checkit} mode; missing otherwise.}
  \item{indoutfilename}{Character: Detailed individual output
  file. Default: "ind.out" in non-\code{checkit} mode; missing
  otherwise.}
  \item{snpoutfilename}{Character: Detailed marker output file. Default:
  "snps.out" in non-\code{checkit} mode; missing otherwise.}
  \item{ethnicfilename}{Character: Output file for average ethnicity at each marker.}
% Extra parameters
  \item{noxdata}{Logical: Indicates presence (TRUE) of X Chromosome data.}
  \item{fakespacing}{Numeric: Spacing between fake markers in
  Morgans. Default: 0.01 (in Morgans).}
  \item{thxpars}{Numeric: Vector of length 3. Sets the initial
    parameters for the prior distribution for theta-X. Default: c(40.0,
    1.0, 10.0).}
  \item{thpars}{Numeric: Vector of length 2. Sets the initial parameters
    for the prior distribution for thea. Default: c(1.0, 5.0).}
  \item{lampars}{Numeric: Vector of length 2. Sets the initial
    parameters for the prior distribution of lambda. Default: c(1.0, 0.1).}
  \item{lamxpars}{Numeric: Vector of length 2. Sets the initial
    parameters for the prior distribution of lambda-X. Default: c(1.0, 0.1).}
  \item{markersim}{Integer: This is the marker of the disease allele, -1
    indicates none. Default: -1.}
  \item{simnumindivs}{Integer: Simulate data with \code{simnumindivs}
    number. Half will be cases, half controls. Half are female, half are
    male. Default: -1.}
  \item{risksim}{Numeric: Risk used to generate data in simulation mode. Default: 1.0.}
  \item{tauscal}{Numeric: Vector of length 2. Initial values of
    tau(African) and tau(European). Default: c(100, 100).}
  \item{wrisk}{Numeric: Allows the model to have weights. Default: 1.0.}
  \item{lrisk}{Numeric: In \code{checkit} mode, leave one marker out in
    turn and this is the risk that we use. Default: -1.0.}
  \item{controlrisk}{Numeric: Risk for controls. Default: 1.0.}
  \item{risk2}{Numeric: Risk for ethnic homozygotes; if NA, \code{risk}\^2 is used.}
  \item{taulsdev}{Numeric: Prior standard deviation for African and
    European tau values. Default: 0.5.}
  \item{taulmean}{Numeric: Prior mean for log10(tau) for both Africans
    and Europeans. Default: 2.0.}
  \item{allmale}{Logical: Data contain only males. Default: FALSE.}
  \item{allcases}{Logical: Data contain only cases. Default: FALSE.}
  \item{usecontrols}{Logical: Include controls in analysis (TRUE). Default: TRUE}
  \item{pbfmodern}{Logical: Publish ancestral allele frequency
    estimates. If TRUE this allows publication of modern allele frequencies. Default: FALSE}
  \item{pbxname}{Character: Debug file for a particular marker. Default: missing.}
  \item{genotoyoutfilename}{Character: Genotype data file generated in
    simulation mode. Default: missing.}
  \item{indtoyoutfilename}{Character: Individual data file generated in
    simulation mode. Default: missing.}
  \item{unknowngender}{Character: Designation for unknown
    gender. Default: "U"}
  \item{ethnic.control}{Character: Name of ethnicity file for
    controls. Not fully implemented.}
  \item{parameter.name}{Character: Name of extra parameters to
    include. Should correspond to \code{value}. See Details.}
  \item{value}{Character Vector: Values of extra parameters to
    include. Should correspond to \code{parameter.name}. See Details.}
  \item{run}{Logical: Run ANCESTRYMAP from R. Will only work on UNIX /
    LINUX systems. Must have ANCESTRYMAP installed and in the PATH.}
  \item{in.prefix}{Character: Source directory for input files.}
  \item{out.prefix}{Character: Destination directory for output files.}
  \item{remove.files}{Logical: Remove temporary files when finished. --Deprecated?}
}

\details{Parameters that are missing (NA) will not be included in the
  parameter file. Parameters that are not included may or may not have a
  default in the ANCESTRYMAP software. Defaults are listed in the
  \code{arguments} section above.

  Additional parameters may be added using the \code{parameter.name} and
  \code{value} arguments. This is the proper way to include future / undocumented /
  obscure parameters at this point.

  In the ANCESTRYMAP software documentation there are 4 different
  example parameter files used for data checking, simulations, and final
  analysis. There are presets available for each of these settings (see
  \code{preset}). Note that in order to add parameters to the default
  list of parameters, you need to include them with the
  \code{parameter.name} and \code{value} arguments.
  
  required parameters, what files are needed/created}

\value{An amap object (see \code{\link{read.amap}})}

\references{See Patterson et al, Methods for High-Density Admixture Mapping of Disease Genes, American Journal of Human Genetics
  (74:5, pp 979-1000)}

\author{Randy Johnson}

\note{ANCESTRYMAP software used by this package was written by Patterson et al. See references.}

\seealso{See also \code{\link{read.amap}} for reading in output files
  from ANCESTRYMAP. Web pages with more information are available at}

\examples{#See amap.pdf}

\keyword{manip}
\keyword{interface}
