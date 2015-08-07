# run.sims.R
# Run the simulated data sets generated in sample.populations.R
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboraotry
# Leidos Biomedical Research, Inc
# Split off of sample.populations.R on May 9, 2014
# Last Modified May 28, 2014


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

    pops <- '/data/johnsonra/mald'
}else{
    library(ALDsuite)

    pops <- 'pops'
}

load(paste(pops, '/anc', seed, '.', OR, '.', n, '.RData', sep = ''))

# make sure we have Pm.priors defined
## if(class(try(load(file = 'Pm.prior.geno.RData'))) == 'try-error')
## {
##     Pm.prior.geno <- setup.prior(hapmap$rs, c('YRI', 'CEU'), cM.linked = 0.1, window = 0.1, maxpcs = NULL,
##                                  dev = FALSE)
##     save(Pm.prior.geno, file = 'Pm.prior.geno.RData')
## }

if(class(try(load(file = paste(pops, '/Pm.prior.haps.RData', sep = '')))) == 'try-error')
{
    Pm.prior.haps <- setup.prior(hapmap$rs, c('YRI', 'CEU'), cM.linked = 0.1, window = 0.1,
                                 phased = TRUE, maxpcs = NULL, dev = FALSE)
    save(Pm.prior.haps, file = paste(pops, '/Pm.prior.haps.RData', sep = ''))
}


######### Run the analysis #########

system.time(out.phased <- admixture(haps = haps, gender = ifelse(anc$male == 1, 'M', 'F'), Pm.prior = Pm.prior.haps,
                             chr = hapmap$chr[hapmap$rs %in% names(Pm.prior.haps)],
                             pos = hapmap$cM[hapmap$rs %in% names(Pm.prior.haps)],
                             indiv.id = anc$pid, burn = 100, iter = 200,
                             marker.id = names(Pm.prior.haps), pop.id = c('YRI', 'CEU'), rand.seed = seed,
                             cores = 1, lambda = 6))

## dads <- 1:(dim(haps)[1] / 2) * 2
## moms <- dads - 1
## geno <- haps[moms,] + haps[dads,]

## system.time(out.unphased <- admixture(haps = NULL, geno = geno, gender = rep('M', dim(haps)[1] / 2),
##                              Pm.prior = Pm.prior.geno, chr = hapmap$chr[hapmap$rs %in% names(Pm.prior.geno)],
##                              pos = hapmap$cM[hapmap$rs %in% names(Pm.prior.geno)],
##                              indiv.id = pick.me, burn = 50, iter = 100,
##                              marker.id = names(Pm.prior.geno), pop.id = c('YRI', 'CEU'), rand.seed = seed,
##                              cores = 1))

save.image(file = paste(paste(pops, '/anc', seed, '.', OR, '.', n, '.RData', sep = '')))


######### How did we do? #########

gammas.sub <- gammas[,dimnames(out.phased$gammas)[[2]],,]
phased.accuracy <- out.phased$gammas - gammas.sub

tmp <- table(as.logical(gammas.sub[,,,1]), round(out.phased$gammas[,,,1])) / length(gammas.sub[,,,1])

if(FALSE)
{
    par(mfrow = c(2, 1))

    # observed
    ALDplot.gammas(out.phased$gammas, indiv = 1)

    # expected
    ALDplot.gammas(gammas.sub, indiv = 1)
}

byMarker <- apply(phased.accuracy[,,,1], 2, mean)

## # unphased
## afri <- with(out.unphased, 2*gammas[,,,1] + gammas[,,,2])
## euro <- with(out.unphased, 2*gammas[,,,3] + gammas[,,,2])

## # phased
## afri <- with(out.phased, gammas[,,1,1] + gammas[,,2,1])
## euro <- with(out.phased, gammas[,,1,2] + gammas[,,2,2])

save.image(file = paste(pops, '/anc', seed, '.', OR, '.', n, '.RData', sep = ''))
