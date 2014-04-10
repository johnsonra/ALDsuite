# samples.R
# Sampling functions
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created September 11, 2012
# Last Modified April 7, 2014


##########################################################################
# Sample gammas for a given orientation (forward/reverse) defined by ord #
##########################################################################

# i = make calculations for the ith individual
# ord = order of markers in increasing (forward) or decreasing (reverse) chromosomal order
# geno = matrix of genotypes in {0, 1, 2} (one row per individual, one column per marker)
# P = matrix of reference allele frequencies (one row per SNP, one column per population)
# lambda = matrix of lambda values (one row per individual, one column per parent)
# d = genetic distances
# chr = chromosome numbers for each marker
# Ajk = array of ancestry proportions (dim1 = individual, dim2 = marker, dim3 = parent, dim4 = population)
# Ak = array of ancestry proportions (dim1 = individual, dim2 = parent, dim3 = population)
# direction = forward or reverse
# returns an array of ancestral probabilities with the following dimesions:
#   i - individual
#   j - marker
#   c - chromosome / parent
#   k - population

sample.G <- function(gammas, lambda, d)
{

    G <- gammas # to get the correct format...
    tmp <- apply(gammas, 1:3, rmultinom, n = 1, size = 1)
    for(k in 1:dim(G)[4])
    {
        G[,,,k] <- tmp[k,,,]
    }

    ### smooth out markers with transitions on both sides to reduce noise ###
    doubles <- array(FALSE, dim(G), dimnames(G))
    doubles[,-1,,] <- G[,-1,,] != G[,-dim(G)[2],,] # ...observed recombination

    # identify double transitions -- two in a row (i.e. ancestries like 1 1 2 1 1)
    doubles[,-dim(G)[2],,] <- doubles[,-1,,] & doubles[,-dim(G)[2],,]
    doubles[,dim(G),,] <- FALSE

    # get rid of double transitions, conditional on gammas
    for(i in 1:dim(doubles)[1])
    {
        for(k in 1:dim(gammas)[3])
        {
            for(j in which(doubles[i,,k,1]))
            {
                # likelihood that we actually observed a double
                pdouble <- gammas[i,j,k,which(G[i,j,k,] == 1)] *
                           prod(P.recombination(TRUE, lambda[i,k], d[j:(j+1)]), na.rm = TRUE)

                # decide if we actually want to keep it (keep with prob pdouble)
                if(!as.logical(rbinom(1, 1, pdouble)))
                    G[i,j,k,] <- 1 - G[i,j,k,]
            }
        }
    }

    return(G)
}

ML.G <- function(gammas)
{
    for(i in 1:dim(gammas)[1])
     {
        for(j in 1:dim(gammas)[2])
        {
            for(k in 1:dim(gammas)[3])
            {
                best <- 1
                ties <- 1
                kept <- 1
                decision.pt <- runif(1)

                for(l in 2:dim(gammas)[4])
                {
                    if(gammas[i,j,k,l] == gammas[i,j,k,best])
                        ties <- ties + 1

                    # if we find a more likely choice, keep it
                    if(gammas[i,j,k,l] > gammas[i,j,k,best])
                    {
                        best <- l
                        ties <- 1
                    }
                }

                if(ties > 1)
                {
                    # keep track of the nth tie
                    n <- 0

                    for(l in 1:dim(gammas)[4])
                    {
                        # if we have ties, pick a uniform random choice from the mix
                        if(gammas[i,j,k,l] == gammas[i,j,k,best])
                        {
                            if(n / ties < decision.pt)
                            {
                                best <- l
                            }
                        }

                        n <- n + 1
                    }
                }

                for(l in 1:dim(gammas)[4])
                {
                    if(l == best)
                    {
                        gammas[i,j,k,l] <- 1
                    }else{
                        gammas[i,j,k,l] <- 0
                    }
                }
            }
        }
    }

    return(gammas)
}

#################
# Sample lambda #
#################

sample.lambda <- function(k, G, Ak, lambda, Pm.prior, distances, chr, alpha, beta, phased, dev)
{
    d.sum <- sum(distances, na.rm = TRUE)

    for(i in 1:dim(G)[1])
    {
        # reset variables for new individual
        n.crossovers <- 0

        begin <- 1 # location of beginning of current segment
        anc <- which(G[i,1,k,] == 1) # ancestry of current segment
        xover <- FALSE

        for(j in 1:dim(G)[2])
        {
            # move along the chromosome until we reach an observed recombination
            if(j == dim(G)[2]) # don't run off the end of the array!
            {
                xover <- TRUE
            }else if(!all(G[i,j,k,] == G[i,j+1,k,]) | chr[j] != chr[j + 1]){
                xover <- TRUE
            }

            if(xover & phased)
            {
                # add the window sizes into the equation (should make minimal difference) ######## fix this!
                if(j > begin)
                {
                    size <- sum(distances[begin:j], na.rm = TRUE) + Pm.prior[[j]]$model[[anc]]$d
                }else{
                    size <- Pm.prior[[j]]$model[[anc]]$d # size of the first window
                }

                # Sample the number of crossovers in this segment
                pval <- runif(1)
                n.sampled <- 0
                lambda.dot <- lambda[i,k]*size / 100
                psum <- ppois(n.sampled, lambda.dot) * Ak[i,k,anc]^n.sampled *
                        exp((1 - Ak[i,k,anc])*lambda.dot)

                while(pval > psum)
                {
                    # see labbook entry from March 13, 2014
                    # titled "Porbability of number of crossovers in a chromosomal segment"
                    n.sampled <- n.sampled + 1
                    psum <- psum + (ppois(n.sampled, lambda.dot) - ppois(n.sampled - 1, lambda.dot)) *
                                   Ak[i,k,anc]^n.sampled * exp((1 - Ak[i,k,anc])*lambda.dot)
                }

                n.crossovers <- n.crossovers + n.sampled

                # don't forget the crossover that resulted in the new segment
                # unless we're moving to a new chromosome
                if(chr[j] == chr[j + 1] & j != dim(G)[2])
                    n.crossovers <- n.crossovers + 1

                # reset variables for new segment
                if(j < dim(G)[2])
                {
                    begin <- j + 1
                    anc <- which(G[i,j + 1,k,] == 1)
                    xover <- FALSE
                }
            }else if(xover & !phased){
                # add the window sizes into the equation below (should make minimal difference)
                if(j > begin)
                {
                    size <- sum(distances[begin:j], na.rm = TRUE) * 2
                }else{
                    size <- 0 # size of the first window will be added in below
                }

                # Sample the number of crossovers in this segment --
                # for each chromosome (may have same or different ancestral states)
                for(anc in as.numeric(c(substr(names(anc), 2, 2), substr(names(anc), 3, 3))))
                {
                    pval <- runif(1)
                    n.sampled <- 0
                    lambda.dot <- lambda[i,k] * size / 100 #Pm.prior[[j]]$model[[anc]]$d
                    psum <- ppois(n.sampled, lambda.dot) * Ak[i,k,anc]^n.sampled *
                                  exp((1 - Ak[i,k,anc])*lambda.dot)

                    while(pval > psum)
                    {
                        # see labbook entry from March 13, 2014
                        # titled "Porbability of number of crossovers in a chromosomal segment"
                        n.sampled <- n.sampled + 1
                        psum <- psum + (ppois(n.sampled, lambda.dot) - ppois(n.sampled - 1, lambda.dot)) *
                                Ak[i,k,anc]^n.sampled * exp((1 - Ak[i,k,anc])*lambda.dot)
                    }

                    n.crossovers <- n.crossovers + n.sampled
                }

                # don't forget the crossover that resulted in the new segment
                # unless we're moving to a new chromosome
                if(chr[j] == chr[j + 1] & j != dim(G)[2])
                    n.crossovers <- n.crossovers + 1

                # reset variables for new segment
                if(j < dim(G)[2]) # (if we aren't at the end of the chrmosome already)
                {
                    begin <- j + 1
                    anc <- which(G[i,j + 1,k,] == 1)
                    xover <- FALSE
                }
            }
        }

        # sample new value for lambda
        lambda[i,k] <- rgamma(1, shape = alpha + n.crossovers, rate = beta + d.sum * (2 - phased))
    }

    return(lambda[,k])
}

## sample.crossovers <- function(G, A, lambda, distances, chr)
## {
##     .Call("sample_crossovers", as.integer(G), as.numeric(A), as.numeric(lambda),
##           as.numeric(distances), as.integer(chr))
## }

## sample.crossoversR <- function(G, A, lambda, distances)
## {

##     n <- dim(G)[1] # number of individuals
##     m <- dim(G)[2] # number of loci
##     L <- dim(G)[3] # number of populations (parents taken out at this stage

##     # "declare" these to avoid an error on the first time through the loop
##     Gcurr <- NULL

##     n.crossovers <- rep(0, n)

##     for(i in 1:n)
##     {
##         for(j in 1:m)
##         {
##             # ancestry at each locus
##             Gprev <- Gcurr

##             Gcurr <- G[i,j,] %*% 1:2

##             # skip over the first marker of each chromosome
##             if(is.na(distances[j]))
##                 next

##             # key probabilities
##             Precomb <- P.recombination(TRUE, lambda[i], distances[j])
##             P0crossovers.marg <- ppois(0, lambda[i] * distances[j])

##             A. <- 0
##             for(l in 1:L)
##             {
##                 ## A. <- apply(A * G[,j,1,], 1, sum)
##                 A. <- A. + A[i,l] * G[i,j,l]
##             }

##             P0crossovers <- P0crossovers.marg / ((1 - Precomb) + Precomb*A.)

##             # if we observed more than 0 crossovers, count them up
##             if(runif(1) > P0crossovers | Gprev != Gcurr)
##             {
##                 # decide if there was a recombination
##                 if(runif(1) <= Precomb | Gprev != Gcurr)
##                 {
##                     q.crossovers <- runif(1) * Precomb # q in (0, Precomb)
##                     odd <- 1 # odd number of crossovers
##                 }else{
##                     q.crossovers <- runif(1) * (1 - Precomb - P0crossovers.marg)
##                     odd <- 0 # even number of crossovers
##                 }

##                 # count up crossovers for different scenarios
##                 quant <- 1 - odd # start at 0 if odd, 1 if even
##                 p.sum <- 0

##                 while(p.sum < q.crossovers & quant < 10000)
##                 {
##                     p.sum <- p.sum + dpois(quant * 2 + odd, lambda[i]*distances[j])
##                     quant <- quant + 1
##                 }

##                 quant <- (quant - 1) * 2 + odd

##                 # add up counts :)

##                 n.crossovers[i] <- n.crossovers[i] + quant
##             }
##         }
##     }

##     return(n.crossovers)
## }

## sample.lambda <- function(G, Ak, lambda, distances, chr, alpha, beta, P.Aeq1, dev)
## {
##     d.sum <- sum(distances, na.rm = TRUE) # total number of centimorgans

##     # number of markers and intervals
##     m <- dim(G)[2]
##     mint <- m - 1

##     # observed recombinations
##     obsRecomb <- !apply(G[,-1,,] == G[,-m,,], 1:3, all)

##     # sampled crossovers
##     n.crossovers <- 0

##     for(j in 1:mint) # working with intervals now...not markers
##     {
##         for(c in 1:2)
##         {
##             # probability of more than one crosover
##             tmpA <- Ak[,c,][as.logical(G[,j,c,])]
##             tmpR <- exp(-lambda[,c] * distances[j+1])
##             p.onePlus <- ifelse(obsRecomb[,j,c], 1, tmpA * (1 - tmpR) / (tmpA * (1 - tmpR) + tmpR))

##             p0 <- ppois(0, lambda[,c] * distances[j+1])
##             q1p <- runif(length(p0), min = p0)

##             samp <- qpois(q1p, lambda[,c] * distances[j+1])

##             n.crossovers <- n.crossovers + ifelse(rbinom(length(p.onePlus), 1, p.onePlus), samp, 0)
##         }
##     }

##     # sample new values for lambda
##     lambda.new <- rgamma(length(n.crossovers), shape = alpha + n.crossovers,
##                          rate = beta + 2*sum(d, na.rm = TRUE)) / 2

##     return(lambda.new)
## }

## # note that this is for each parent's contributed chromosome -- min has to be at least 1
## # ...or discretely 0 if the individual is unadmixed
## sample.lambda.k <- function(G0, G, A, distances, chr, alpha, beta, P.Aeq1, dev)
## {
##     # make sure we don't have any 0's or 1's
##     bad <- apply(A, 1, prod) == 0
##     if(any(bad))
##         A[bad,] <- A[bad,] - A[bad,] * 0.01 + (1 - A[bad,]) * 0.01

##     d.sum <- sum(distances, na.rm = TRUE) # total number of centimorgans

##     m <- dim(G)[2] # get number of loci (there will be m-1 intervals)

##     # make sure we have an over-all change in ancestry before we call this a crossover!
##     change <- matrix(FALSE, nrow = dim(G0)[1], ncol = m - 1)

##     for(l in 1:dim(G0)[3])
##     {
##         change <- change | (G0[,-m,l] != G0[,-1,l])
##     }

##     # count departures from this state (don't include transitions between chromosomes
##     transitions <- array(0, dim = dim(A), dimnames = dimnames(A))

##     for(l in 1:dim(G)[3])
##     {
##         transitions[,l] <- (G[,-m,l] != G[,-1,l] & change) %*% (chr[-m] == chr[-1])
##     }

##     # use number of observed transitions to infer # crossovers
##     n.crossovers <- apply(transitions / (1 - A), 1, sum)

##     # put a lower bound on this, based on A
##     tmp <- apply(A, 1, min) # the minimum will be according to the smallest value in A
##     min.lambda <- (log(1) - log(tmp)) / log(2)

##     # update lambda conditional on posterior distribution of lambda and min.lambda
##     lambda.new <- rep(NA, dim(A)[1])

##     for(i in 1:dim(A)[1])
##     {
##         if(rbinom(1, 1, max(P.Aeq1[i,])))
##         {
##             # no admixture inferred in this sample? '0' is bad...just don't set a minimum?
##             lambda.new[i] <- rgamma(1, shape = alpha + n.crossovers[i], rate = beta + d.sum)
##         }else{
##             p.min <- pgamma(min.lambda[i], shape = alpha + n.crossovers[i], rate = beta + d.sum)
##             p.samp <- runif(1, min = p.min, max = 1)

##             lambda.new[i] <- qgamma(p.samp, shape = alpha + n.crossovers[i], rate = beta + d.sum)
##         }
##     }

##     return(lambda.new)
## }
