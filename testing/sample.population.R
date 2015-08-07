# sample.population.R
# Sample a population for analysis (from make.population.R)
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Split from make.population.R November 26, 2013
# Last Modified May 9, 2014


# add "--args <seed> <OR> <n>" to the end of the command
args <- as.numeric(commandArgs(TRUE))
if(length(args) != 3)
    stop('Wrong number of arguments')

seed <- args[1]
OR <- args[2]
n <- args[3]

if(!require(ALDdata))
{
    library(hwde, lib.loc = '~/lib')
    library(ALDdata, lib.loc = '~/lib')
    library(ALDsuite, lib.loc = '~/lib')
    library(rhdf5, lib.loc = '~/lib')

    pops <- '/data/johnsonra/mald'

}else{
    library(ALDsuite)
    library(rhdf5)

    pops <- 'pops'
}

# set seed
set.seed(seed)

# Load prior data
load(paste(pops, 'anc.RData', sep = '/'))
load(paste(pops, 'hapmap.RData', sep = '/'))

# pick anchor markers
if(class(try(load(file = paste(pops, 'Pm.prior.haps.RData', sep = '/')))) == 'try-error')
{
    Pm.prior.haps <- setup.prior.dev(hapmap$rs, c('YRI', 'CEU'), cM.linked = 0.05, window = 0.05,
                                 phased = TRUE, maxpcs = NULL)
    save(Pm.prior.haps, file = paste(pops, 'Pm.prior.haps.RData', sep = '/'))
}

# pick a risk allele (we've chosen these to be aims, so it should be a decent one)
risk.allele <- sample(names(Pm.prior.haps), size = 1)
risk.allele.j <- which(hapmap$rs == risk.allele)

# get alleles from HDF5 files
files <- system(paste('ls', pops, '| grep .h5'), intern = TRUE)

groups <- gsub('.h5', '', gsub('Population', '', files), fixed = TRUE)

files <- paste(pops, '/', files, sep = '')

for(i in 1:length(groups))
{
    tmp <- h5read(files[i], paste(groups[i], '/gammas', sep = ''), index = list(NULL, risk.allele.j, NULL, 1))
    tmp <- apply(tmp, 1, sum)

    anc$anc[anc$group == i] <- tmp
}

######### set up analysis #########
# don't sample by allele, sample by ancestry!!!!! -- model isn't quite right either...
anc$risk <- OR^anc$anc

# pick cases
cases <- sample(1:dim(anc)[1], size = n, prob = anc$risk)

# pick controls
conts <- sample((1:dim(anc)[1])[-cases], size = n)

# sampled subset
anc$case <- FALSE
anc$case[cases] <- TRUE
pick.me <- c(cases, conts)

anc <- anc[pick.me,]

alleles <- array(NA, dim = c(n*2, dim(hapmap)[1], 2))
gammas <- array(NA, dim = c(n*2, dim(hapmap)[1], 2, 2))
for(i in 1:length(groups))
{
    if(any(anc$group == i))
    {
        alleles[which(anc$group == i),,] <- with(subset(anc, group == i),
                          h5read(files[i], paste(groups[i], '/alleles', sep = ''),
                                 index = list(id, NULL, NULL)))

        gammas[which(anc$group == i),,,] <- with(subset(anc, group == i),
                         h5read(files[i], paste(groups[i], '/gammas', sep = ''),
                                index = list(id, NULL, NULL, NULL)))
    }
}

# create haps variable
dads <- 1:(dim(alleles)[1]) * 2
moms <- dads - 1

haps <- matrix(nrow = dim(alleles)[1] * 2, ncol = dim(alleles)[2],
               dimnames = list(rep(dimnames(alleles)[[1]], each = 2), dimnames(alleles)[[2]]))
haps[moms,] <- alleles[,,1]
haps[dads,] <- alleles[,,2]


############
# Dimnames #
############

anc$pid <- paste(anc$group, anc$id, sep = '.')

dimnames(gammas) <- list(anc$pid, hapmap$rs, c('Mother', 'Father'), c('YRI', 'CEU'))

dimnames(alleles) <- list(anc$pid, hapmap$rs, c('Mother', 'Father'))
dimnames(haps) <- list(rep(anc$pid, each = 2), hapmap$rs)


########
# Save #
########

save(haps, gammas, Pm.prior.haps, n, OR, seed, anc, hapmap, alleles,
     file = paste(pops, '/anc', seed, '.', OR, '.', n, '.RData', sep = ''))
