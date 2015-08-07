# asw.R
# Run ASW data through ALDsuite
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Copied from asw.R on October 18, 2013
# Last Modified October 31, 2013

library(ALDsuite)
library(ALDdata)
source('../R/helpers.R')

if(class(try(load('mex.RData'))) == 'try-error')
{
    nMarkers <- 1000
    chr <- 22

    # get data for one chromosome
    mex <- get.phased(chr, 'MEX', clean.duos = TRUE)

    # subsample
    set.seed(2380947)
    mex <- mex[sample(1:dim(mex)[1], nMarkers,),]

    rownames(mex) <- mex$rsID
    mex$rsID <- NULL
    mex$position_b36 <- NULL

    # get reference information from YRI and convert to appropriate format
    data(hapmap)
    hapmap <- subset(hapmap, rs %in% rownames(mex))
    hapmap <- hapmap[order(hapmap$pos),]
    mex <- mex[hapmap$rs,] # make sure we are ordered the same!!!

    geno <- t(categorize(as.matrix(mex), ref = hapmap$ref))

    # generate our PCR prior
    Pm.prior <- setup.prior(hapmap$rs, c('YRI', 'CEU', 'CHB+JPT'), maxpcs = 1)

    # set up the remaining variables for admixture
    gender <- rep('M', dim(geno)[1])
    chr <- rep(chr, dim(geno)[2])
    pos <- hapmap$cM
    build <- 36
    d <- NULL
    burn <- 100
    iter <- 200
    every <- 2
    indiv.id <- rownames(geno)
    marker.id <- colnames(geno)
    pop.id <- c('YRI', 'CEU', 'CHB+JPT')
    lambda <- 6
    tau <- 300
    omega <- NULL
    rand.seed <- 82430
    cores <- 1
    cl <- NULL
    dev <- FALSE
    verbose <- TRUE
    debug <- FALSE
    sex.chr <- 23
    indiv.geno.check <- 0.98
    marker.geno.check <- 0.98
    hwe.thresh <- 1e-4
    bad.indiv <- NULL
    bad.marker <- NULL
    male.het.X <- 1
    fast <- FALSE

    save(geno, Pm.prior, gender, chr, pos, build, d, burn, iter, every, indiv.id, marker.id,
         pop.id, lambda, tau, omega, rand.seed, cores, cl, dev, verbose, debug, sex.chr,
         indiv.geno.check, marker.geno.check, hwe.thresh, bad.indiv, bad.marker, fast,
         file = 'mex.RData')

    rm(list = ls())

    load('mex.RData')
}

system.time(out <- admixture(geno, gender, Pm.prior, chr, pos, burn, iter, every, indiv.id, marker.id,
                             pop.id, lambda, tau, omega, rand.seed, cores, cl, dev, verbose, debug,
                             sex.chr, indiv.geno.check, marker.geno.check, hwe.thresh, bad.indiv,
                             bad.marker, male.het.X, fast))

save(out, file = 'mex_out.RData')
