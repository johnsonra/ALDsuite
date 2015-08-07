# asw.R
# Run ASW data through ALDsuite
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created September 3, 2013
# Last Modified April 4, 2014

library(ALDdata)
library(ALDsuite)

source('../../ALDdata/R/helpers.R')


if(class(try(load('asw0.4.0.RData'))) == 'try-error')
{
    nMarkers <- 20000
    chr2sample <- 20

    # get data for one chromosome
    asw <- get.phased(chr2sample, 'ASW')

    rownames(asw) <- asw$rsID
    asw$rsID <- NULL
    asw$position_b36 <- NULL

    # get reference information from YRI and convert to appropriate format
    data(hapmap)
    hapmap <- subset(hapmap, rs %in% rownames(asw))
    hapmap <- hapmap[order(hapmap$pos),]
    asw <- asw[hapmap$rs,] # make sure we are ordered the same!!!

    geno <- t(categorize(as.matrix(asw), ref = hapmap$ref))
    haps <- t(categorize(as.matrix(asw), ref = hapmap$ref, collapse = FALSE))

    # drop monomorphic markers
    mono <- apply(geno, 2, sum) == 0

    geno <- geno[,!mono]
    haps <- haps[,!mono]
    hapmap <- subset(hapmap, rs %in% colnames(geno))

    # generate our PCR prior on a subsample of the rs numbers
    set.seed(2380947)
    Pm.prior <- setup.prior.dev(hapmap$rs, c('YRI', 'CEU'), maxpcs = NULL, cM.linked = 0.1, window = 0.1,
                                phased = TRUE)
    Pm.prior.geno <- setup.prior.dev(names(Pm.prior), c('YRI', 'CEU'), cM.linked = 0.1, window = 0.1)

    ## system('say done')

    # set up the remaining variables for admixture
    gender <- rep('M', dim(geno)[1])
    chr <- rep(chr2sample, length(Pm.prior))
    pos <- hapmap$cM[hapmap$rs %in% names(Pm.prior)]
    build <- 36
    d <- NULL
    indiv.id <- rownames(geno)
    marker.id <- names(Pm.prior)
    pop.id <- c('YRI', 'CEU')
    lambda <- 6
    tau <- 300
    omega <- NULL
    rand.seed <- 82430
    cores <- 1
    cl <- NULL
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
    nPCs <- 1

    dev <- FALSE
    burn <- 10
    iter <- 20
    every <- 5

    dads <- 1:(dim(haps)[1] / 2) * 2
    moms <- dads - 1
    geno <- haps[moms,] + haps[dads,]

    save.image('asw0.4.0.RData')
}

if(TRUE)
{
    system.time(out.phased <- admixture(Pm.prior = Pm.prior, haps = haps, geno = NULL, gender = gender,
                                        chr = chr, pos = pos, burn = burn, iter = iter, every = every,
                                        indiv.id = indiv.id, marker.id = marker.id, pop.id = pop.id,
                                        lambda = lambda, tau = tau, omega = omega, rand.seed = rand.seed,
                                        cores = cores, cl = cl, dev = dev, verbose = verbose, debug = debug,
                                        sex.chr = sex.chr, indiv.geno.check = indiv.geno.check,
                                        marker.geno.check = marker.geno.check, hwe.thresh = hwe.thresh,
                                        bad.indiv = bad.indiv, bad.marker = bad.marker,
                                        male.het.X = male.het.X, fast = fast))
    save.image(file = 'asw0.4.0.RData')

    system.time(out.unphased <- admixture(Pm.prior = Pm.prior.geno, haps = NULL, geno = geno, gender = gender,
                                          chr = chr, pos = pos, burn = 100, iter = 200, every = every,
                                          indiv.id = indiv.id, marker.id = marker.id, pop.id = pop.id,
                                          lambda = lambda, tau = tau, omega = omega, rand.seed = rand.seed,
                                          cores = cores, cl = cl, dev = dev, verbose = verbose, debug = debug,
                                          sex.chr = sex.chr, indiv.geno.check = indiv.geno.check,
                                          marker.geno.check = marker.geno.check, hwe.thresh = hwe.thresh,
                                          bad.indiv = bad.indiv, bad.marker = bad.marker,
                                          male.het.X = male.het.X, fast = fast))
    save.image(file = 'asw0.4.0.RData')
}

if(FALSE)
{
    colors <- c(rgb(255, 102, 102, maxColorValue = 255), rgb(0, 128, 128, maxColorValue = 255))
    ## colors <- c('blue', 'green3')

    pdf('phased.pdf')
    par(mfrow = c(8, 1), mar = c(.5,2,.5,1))
    for(i in 1:dim(out.phased$gammas)[1])
    {
        afri <- with(out.phased, cbind(gammas[i,,1,1], gammas[i,,2,1]))
        euro <- with(out.phased, cbind(gammas[i,,1,2], gammas[i,,2,2]))

        barplot(rbind(euro[,2], afri[,2], euro[,1], afri[,1]), space = 0, axisnames = FALSE, axes = FALSE,
                col = rep(colors, 2), border = NA)

        lines(c(1, dim(out.phased$gammas)[2]), rep(1, 2))
    }
    dev.off()

    pdf('unphased.pdf')
    par(mfrow = c(8, 1), mar = c(.5,2,.5,1))
    for(i in 1:dim(out.unphased$gammas)[1])
    {
        afri <- out.unphased$gammas[i,,,1]*2 + out.unphased$gammas[i,,,2]
        euro <- 2 - afri

        barplot(rbind(euro, afri), space = 0, axisnames = FALSE, axes = FALSE,
                col = colors, border = NA)
    }
    dev.off()
}


######### Other analyses #########
if(FALSE)
{
    ##### PCAdmix #####

    # beagle files
    bhaps <- cbind('M', colnames(haps), t(haps))
    write.table(bhaps, file = 'admixed.txt', sep = '\t', row.names = FALSE, col.names = FALSE,
                quote = FALSE)

    data(ceu20)
    bceu <- phased[,colnames(haps)]
    bceu <- cbind('M', colnames(bceu), t(bceu))
    write.table(bceu, file = 'ceu.txt', sep = '\t', row.names = FALSE, col.names = FALSE,
                quote = FALSE)

    data(yri20)
    byri <- phased[,colnames(haps)]
    byri <- cbind('M', colnames(byri), t(byri))
    write.table(byri, file = 'yri.txt', sep = '\t', row.names = FALSE, col.names = FALSE,
                quote = FALSE)

    ## system('./PCAdmix3_macOSX -anc yri.txt ceu.txt -adm admixed.txt')

    raw <- read.table('ancestry_output.fbk.txt')

    pcadmix <- array(0, dim = c(dim(out.phased$gammas)[1], dim(raw)[2] - 2, 2, 2),
                     dimnames = list(dimnames(out.phased$gammas)[[1]], NULL, c('Mother', 'Father'),
                                     c('YRI', 'CEU')))

    i <- 1:dim(out.phased$gammas)[1] * 4

    for(j in 3:dim(raw)[2])
    {
        pcadmix[,j-2,1,1] <- raw[i-3,j]
        pcadmix[,j-2,1,2] <- raw[i-2,j]
        pcadmix[,j-2,2,1] <- raw[i-1,j]
        pcadmix[,j-2,2,2] <- raw[i  ,j]
    }

    # calculate widths
    widths <- rep(0, dim(pcadmix)[2])
    windows <- read.table('ancestry_output.markers.txt', sep = '\t', stringsAsFactors = FALSE)$V2

    for(j in 1:length(widths))
    {
        snps <- strsplit(windows[j], ' ', fixed = TRUE)[[1]]
        tmp <- range(subset(hapmap, rs %in% snps)$cM)
        widths[j] <- tmp[2] - tmp[1]
    }

    widths <- widths + 0.01

    pdf('PCAdmix.pdf')
    par(mfrow = c(8, 1), mar = c(.5,2,.5,1))
    for(i in 1:dim(pcadmix)[1])
    {
        afri <- cbind(pcadmix[i,,1,1], pcadmix[i,,2,1])
        euro <- cbind(pcadmix[i,,1,2], pcadmix[i,,2,2])

        barplot(rbind(euro[,2], afri[,2], euro[,1], afri[,1]), space = 0, axisnames = FALSE, axes = FALSE,
                col = rep(colors, 2), border = NA, width = widths)

        ## lines(c(0, length(widths)), rep(1, 2))
        lines(c(0, sum(widths)), rep(1, 2))
    }
    dev.off()

    ##### MULTIMIX #####

    # waco formatting for haplotypes...I'm sure this helps it run a lot faster though...
    waco <- function(x) ifelse(x == 1, "1 0 0", "0 0 1")

    mmixhaps <- cbind(hapmap$rs, hapmap$pos, hapmap$ref, hapmap$var, apply(haps, 1, waco))
    write.table(mmixhaps, file = 'multimix/data/Samples/haplos/haplos_chr20', sep = ' ', quote = FALSE,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    data(ceu20)
    write.table(t(phased[,hapmap$rs]), file = 'multimix/data/haplotypes/CEU/chr20.haps', sep = ' ',
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    data(yri20)
    write.table(t(phased[,hapmap$rs]), file = 'multimix/data/haplotypes/YRI/chr20.haps', sep = ' ',
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    legnd <- cbind(hapmap$rs, hapmap$pos, hapmap$ref, hapmap$var)
    colnames(legnd) <- c('rs', 'position', 'a0', 'a1')
    write.table(legnd, file = 'multimix/data/haplotypes/legend_files/chr20.legend', sep = ' ',
                row.names = FALSE, quote = FALSE)

    n <- dim(hapmap)[1]
    map <- cbind(hapmap$pos,
                 with(hapmap, sapply(c(7.77, (cM[-1] - cM[-n]) / (pos[-1] - pos[-n]) * 1e6), min, 300)),
                 hapmap$cM - min(hapmap$cM))
    colnames(map) <- c('position', 'COMBINED_rate(cM/Mb)', 'Genetic_Map(cM)')
    write.table(map, file = 'multimix/data/haplotypes/genetic_maps/chr20.map', sep = ' ',
                row.names = FALSE, quote = FALSE)

    #./MULTIMIX_MCMC -dir ./data/ -o ./output/ -sourcepops CEU/ YRI/ -misfit 0.95 0.05 0.05 0.95 -chr.list 20

    # set window widths
    snps.per.window <- read.table('multimix/output/snps_per_window.txt')$V3
    widths <- rep(0, length(snps.per.window))

    at <- 1

    for(j in 1:length(widths))
    {
        tmp <- range(hapmap[at:(at + snps.per.window[j] - 1),]$cM)
        widths[j] <- tmp[2] - tmp[1]

        at <- at + snps.per.window[j]
    }

    widths <- widths + 0.01

    pdf('MULTIMIX.pdf')
    par(mfrow = c(8, 1), mar = c(.5,2,.5,1))
    for(i in 1:(dim(haps)[1] / 2))
    {
        hap1 <- read.table(paste('multimix/output/Local_samples_indiv', i, '_haplo1.txt', sep = ''))
        hap2 <- read.table(paste('multimix/output/Local_samples_indiv', i, '_haplo2.txt', sep = ''))

        afri <- cbind(hap1$V4, hap2$V4) / 100
        euro <- cbind(hap1$V3, hap2$V3) / 100

        barplot(rbind(euro[,2], afri[,2], euro[,1], afri[,1]), space = 0, axisnames = FALSE, axes = FALSE,
                col = rep(colors, 2), border = NA, width = widths)

        ## lines(c(0, length(widths)), rep(1, 2))
        lines(c(0, sum(widths)), rep(1, 2))
    }
    dev.off()
}
