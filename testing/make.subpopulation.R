# make.subpopulation.R
# Randall Johnson
# Taken from make.population.R to help with memory issues
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Taken from make.population.R May 5, 2014
# Last Modified May 6, 2014


## n = how many we want
## N = how many to simulate (to pick from)
## lastn = last identifier (to be added to 1:n for anc$id)

make.subpopulation <- function(n, N, lastn, ceu, yri, hapmap)
{
    hapmap <- subset(hapmap, rs %in% colnames(ceu))
    nMarkers <- dim(hapmap)[1]

    ###############################################
    # Assign values for each individual's parents #
    ###############################################

    # Save matricies of these to speed things up
    # I'm using the transpose of the matrix because it will write more quickly this way (in my loop below)
    haps <- matrix(ncol = 2*N, nrow = nMarkers, dimnames = list(hapmap$rs))
    gammas <- matrix(ncol = 2*N, nrow = nMarkers, dimnames = list(hapmap$rs))

    for(i in 1:(2*N))
    {
        tmp <- 0
        at <- 1

        while(at <= nMarkers)
        {
            tmp <- tmp + 1
            G <- sample(1:2, 1, prob = c(0.82, 0.18))

            if(G == 1)
            {
                ancHap <- yri[sample(1:dim(yri)[1], 1),]
            }else{
                ancHap <- ceu[sample(1:dim(ceu)[1], 1),]
            }

            hop <- rexp(1, 6) * 100
            to <- with(hapmap, max(c((1:length(cM))[cM < cM[at] + hop], at))) # end of this block

            gammas[at:to, i] <- G
            haps[at:to, i] <- ancHap[at:to]

            at <- to + 1
        }

        if(any(is.na(gammas[,i]) | is.na(haps[,i])))
        {
            save.image(file = paste('anc', seed, '_error.RData', sep = ''))
            stop()
        }
    }

    ### sample haps to represent the distribution I want ###
    # the results above are too broad, but this might be a little restraining...?

    # observed ancestry
    A <- apply(2 - gammas, 2, mean)
    pA <- dbeta(A, 8.2*7, 1.8*7)
    recombs <- apply(gammas, 2, function(x) sum(x[-1] != x[-length(x)]))

    # pick individuals that I want
    pick.me <- sample(1:(2*N), size = 2*n, prob = pA)

    ### calculate metrics of sampled chromosomes ###

    # indexing haps to only get mom's or dad's chromosome
    moms <- (1:n - 1) * 2 + 1
    dads <- (1:n) * 2

    anc <- data.frame(id = lastn + 1:n,
                      p1.id = paste('p1', 1:n, sep = '-'),
                      p2.id = paste('p2', 1:n, sep = '-'),
                      A1 = A[pick.me][moms],
                      A2 = A[pick.me][dads],
                      gender = sample(c('M', 'F'), size = n, replace = TRUE))
    anc$A0 <- (anc$A1 + anc$A2) / 2

    ### subset gammas and haps -- reformat ###

    haps <- t(haps[,pick.me])
    tmp <- t(gammas[,pick.me])

    gammas <- array(0, dim = c(n, nMarkers, 2, 2),
                    dimnames = list(anc$id, hapmap$rs, c("Mother", "Father"), c("YRI", "CEU")))

    gammas[,,1,1][tmp[moms,] == 1] <- 1
    gammas[,,2,1][tmp[dads,] == 1] <- 1
    gammas[,,1,2][tmp[moms,] == 2] <- 1
    gammas[,,2,2][tmp[dads,] == 2] <- 1

    alleles <- array(0, dim = c(n, nMarkers, 2),
                     dimnames = list(anc$id, hapmap$rs, c("Mother", "Father")))

    alleles[,,1] <- haps[moms,]
    alleles[,,2] <- haps[dads,]


    return(list(anc = anc, gammas = gammas, alleles = alleles))
}


setup.hdf5File <- function()
{
    # create file
    h5createFile(maldPop)

    # create group
    h5createGroup(maldPop, groupName)

    # create anc
    h5createDataset(maldPop, '/anc', c(subn, 6), storage.mode = 'double', chunk = c(subn, 1))

    h5write(1:subn, maldPop, '/anc', start = c(1,1), count = c(subn,1))
    h5write(rep(group, subn), maldPop, '/anc', start = c(1,2), count = c(subn, 1))

    # create gammas and alleles
    h5createDataset(maldPop, paste('/', groupName, '/gammas', sep = ''),
                    c(subn, length(tokeep), 2, 2), storage.mode = 'integer',
                    chunk = c(1, length(tokeep), 2, 2))
    h5createDataset(maldPop, paste('/', groupName, '/alleles', sep = ''),
                    c(subn, length(tokeep), 2), storage.mode = 'integer',
                    chunk = c(1, length(tokeep), 2))
}
