# make.population.R
# Simulate a population from which to sample
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created August 1, 2013
# Last Modified May 8, 2014


##### Globals #####

# sample 1 million individuals -- effective population
Neff <- 1e6

# chromosome
chr <- 20

# number of groups to break Neff individuals into
nGroups <- 100

# number of individuals to sample in each iteration of this run
n <- 10

# number of simulated indivdiuals to create for each iteration of this run
N <- n * 100

# number of individuals to generate in this subpopulation
subn <- Neff / nGroups

# random seed - will be reset conditional upon group later on
set.seed(29834729)

##### Setup #####

# add "--args <group>" to the end of the command
# <group> defines which group we are creating in this run

args <- as.numeric(commandArgs(TRUE))
if(length(args) != 1)
    stop('Wrong number of arguments')

group <- args[1]
groupName <- paste(paste(rep(0, 2 - floor(log10(group))), collapse = ''), group, sep = '') # for 001:999

if(!require(ALDdata))
{
    library(hwde, lib.loc = '~/lib')
    library(ALDdata, lib.loc = '~/lib')
    library(ALDsuite, lib.loc = '~/lib')
    library(rhdf5, lib.loc = '~/lib')
    setwd('~/mald/')

    maldPop <- paste('/data/johnsonra/mald/Population', groupName, '.h5', sep = '')
}else{
    library(ALDsuite)
    library(rhdf5)

    maldPop <- paste('pops/Population', groupName, '.h5', sep = '')
}

source('make.subpopulation.R')

##### Set up Hapmap data #####

data(ceu20)
ceu <- phased

data(yri20)
yri <- phased

data(hapmap)

# pick markers (keep all to start with)
tokeep <- 1:dim(ceu)[2]

# make sure we don't have any monomorphic markers
tmp <- which(apply(ceu[,tokeep], 2, sum) == 0 & apply(yri[,tokeep], 2, sum) == 0)
tokeep <- tokeep[-tmp]

# subset datasets
ceu <- ceu[,tokeep]
yri <- yri[,tokeep]

hapmap <- subset(hapmap, rs %in% colnames(ceu))
hapmap <- hapmap[order(hapmap$pos),]

if(FALSE)
    save(hapmap, file = 'pops/hapmap.RData')

###### create netCDF file for housing our large population ######
setup.hdf5File()

### Set seed for this group ###
seeds <- unique(round(runif(nGroups*2) * 100000))[1:nGroups]

set.seed(seeds[group])

### make the subpopulations of interest ###
lastn <- 0

while(lastn < subn)
{
    tmp <- make.subpopulation(n, N, lastn, ceu, yri, hapmap)

    h5write(with(tmp$anc, cbind(A1, A2, A0, as.numeric(gender == 'M'))), file = maldPop, name = '/anc',
            start = c(lastn + 1, 3), count = c(n, 4))

    h5write(tmp$gammas, file = maldPop, name = paste('/', groupName, '/gammas', sep = ''),
            start = c(lastn + 1, 1, 1, 1), count = c(n, length(tokeep), 2, 2))

    h5write(tmp$alleles, file = maldPop, name = paste('/', groupName, '/alleles', sep = ''),
            start = c(lastn + 1, 1, 1), count = c(n, length(tokeep), 2))

    lastn <- lastn + n
}
