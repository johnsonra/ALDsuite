# likelihoods.R
# Distribution of ancestral state (gamma) given necessary parameters
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created September 1, 2012
# Last Modified February 6, 2014

#########
# Notes #
#########

# - These probabilities are for the HMM step (randomly selecting gammas),
#   so many random variables are treated as fixed and will be updated elsewhere.

# - These probabilities are for *one* chromosome.
#   Each chromosome is treated separately.

# - Error checking on these functions is rather limited, as they are
#   not intended to be used by the end user.

######################
# function variables #
######################

# x = value for which to calculate the probability
# lambda ~= average number of generations since admixture in a parent
# lambda.alt ~= average number of generations since admixture in the *other* parent
# distance = Morgans between previous and current loci
# gamma = ancestry {1,2,3,...} at the current locus
# gamma.alt = ancestry at the current locus for the *other* chromosome
# gamma.prev = value for gamma at previous locus
# gamma.alt.prev = value for gamma at previous locus for the *other* chromosome
# Ai = Individual global ancestry
# Ak = Parent's global ancestry
# Ak.alt = *Other* parent's global ancestry
# P = vector of reference allele frequencies in each ancestral population
# a = Individual genotype at current locus
# alpha = shape parameter for prior on lambda
# beta = rate parameter for prior on lambda
# dir = dirichlet parameter vector for prior on Aks

############################################
# P(recombination == x ; lambda, distance) #
############################################

# x must be of length 1
# lambda and distance may be longer in length (either length of 1 or equal in length)
P.recombination <- function(x, lambda, distance)
{
    .Call("P_recombination", as.integer(x), as.numeric(lambda), as.numeric(distance) / 100)
}

####################################################
# P(gamma == x ; lambda, distance, Ak, gamma.prev) #
####################################################

# x in {1,2,...} (ancestral allele population code)
P.gamma <- function(x, lambda, distance, Ak, gamma.prev = rep(NA, dim(Ak)[2]), p.recomb.T, p.recomb.F)
{
  if(is.null(dim(gamma.prev)))
      gamma.prev <- t(gamma.prev)

  if(is.null(dim(Ak)))
      Ak <- t(Ak)

  # for the first locus in the chain this is simply a function of global ancestry
  if(length(gamma.prev) > 1) # I don't think this is correct...should rely on dimensions, no?
      return(ifelse(is.na(gamma.prev[,x]), Ak[,x],
                    p.recomb.F * gamma.prev[,x] +  # function of recombination
                    p.recomb.T * Ak[,x]))

  return(ifelse(is.na(gamma.prev[x]), Ak[x],
                    p.recomb.F * gamma.prev[x] +  # function of recombination
                    p.recomb.T * Ak[x]))
}

####################################
# P(a == x | gamma, gamma.alt ; P) #
####################################

# x in {0,1,2} (0, 1 or 2 variant alleles)
P.a.gamma.gamma <- function(x, gamma, gamma.alt, P)
{
    if(is.null(dim(P))) # when a vector of length equal to the number of populations
    {                   # don't forget 0-based indexing ___________________________________________
                        #                                                \                         \
                        #                                                 V                         V
       tmp <- .Call("P_a_gamma_gamma1", as.integer(0:2), as.integer(gamma - 1), as.integer(gamma.alt - 1),
                    as.numeric(P))

       return(ifelse(is.na(x), 1, tmp[x + 1]))
   }

    # when an array with dim1 - individual, dim2 - chromosome (mother/father), dim3 - population
    .Call("P_a_gamma_gamma2", as.integer(x), as.integer(gamma - 1), as.integer(gamma.alt - 1),
          as.numeric(P), as.integer(dim(P)))
}

###############################################################################################
# P(gamma == x | a ; P, lambda, lambda.alt, distance, Ak, Ak.alt, gamma.prev, gamma.alt.prev) #
###############################################################################################

P.gammas <- function(geno, P, Pm.prior, lambda, lambdaX, d, chr, A0, Ak, AX, gender, sex.chr, haps, dev)
{
    if(haps)
    {
        combos <- dimnames(P)[[2]]

        gammas <- array(0, c(dim(geno)[1], length(Pm.prior), 2, length(combos)),
                        dimnames = list(dimnames(geno)[[1]], dimnames(P)[[1]],
                                        c('Mother', 'Father'), combos))
    }else{
        combos <- character()

        for(k1 in 1:dim(P)[2])
            for(k2 in k1:dim(P)[2])
                combos <- c(combos, paste('g', k1, k2, sep = ''))

        gammas <- array(0, c(dim(geno)[1], length(Pm.prior), 1, length(combos)),
                        dimnames = list(dimnames(geno)[[1]], dimnames(P)[[1]],
                                        'Empty', combos))
    }


    for(chr.curr in unique(chr))
    {
        ## if(haps)
        ## {
        ##     # probability of g given previous state
        ##     pag <- array(0, c(dim(geno)[1], sum(chr == chr.curr), 2, dim(P)[2]),
        ##                  dimnames = list(dimnames(geno)[[1]], names(Pm.prior)[chr == chr.curr],
        ##                                  c('Mother', 'Father'), dimnames(P)[[2]]))

        ##     # use these for calculations below (saves time)
        ##     pr <- pag[,,,1]
        ##     pgNeqg <- pag
        ## }

        # use for sex chromosome calculations
        if(chr.curr == sex.chr)
        {
            combosX <-  which(substr(combos, 2, 2) == substr(combos, 3, 3))
            M <- which(gender == 'M')
            het <- which(gender == 'F') # which individuals can be hets
        }else{
            het <- 1:dim(geno)[1] # anyone can be a het on the autosomes
        }

        ord <- (1:length(chr))[chr == chr.curr]

        ### Probabilities of recombination and of a | gamma ###
        # Forward chain #
        for(j in ord)
        {
            # this is the chromosome specific index (j is relative to the entire genome)
            jj <- j - min(ord) + 1

            ### P(gammas | a) ###
            if(haps)
            {
                gammas[,j,,] <- P.a.gamma(gammas[,max(j-1, 1),,], geno, Pm.prior[[j]]$model,
                                          P[j,], Ak, lambda, d[j])
            }else{
                # -- gammas.prev is ignored on the first time through --
                gammas[,j,1,] <- P.gammas.a(geno[,j], geno, Pm.prior[[j]], P[j,], A0, Ak,
                                            lambda, d[j], gammas[,max(j - 1, 1),1,])
            }
        }

        g.prev <- gammas[,1,,] # ignored, but this is the correct format...

        # Reverse Chain #
        for(j in ord[length(ord):1])
        {
            # this is the chromosome specific index (j is relative to the entire genome)
            jj <- j - min(ord) + 1

            if(haps)
            {
                g.prev <- P.a.gamma(g.prev, geno, Pm.prior[[j]]$model, P[j,], Ak, lambda, d[j+1])
                gammas[,j,,] <- gammas[,j,,] * g.prev
            }else{
                if(jj == length(ord))
                {
                    g.prev <- gammas[,1,1,, drop = FALSE]
                }

                g.prev[,1,1,] <- P.gammas.a(geno[,j], geno, Pm.prior[[j]], P[j,], A0, Ak,
                                          lambda, d[j+1], g.prev[,1,1,])
            }
        }

        # normalize the two
        sum2one <- function(x) x / sum(x)

        # i = 3, j = 454
        tmp <- apply(gammas, 1:3, sum2one)
        gammas[,,,1] <- tmp[1,,,]
        gammas[,,,2] <- tmp[2,,,]
    }

    return(gammas)
}


# assuming a in {0,1} because this should only be run on X chromosome markers for males
#   --> only one chromosome to be working with here!
P.gamma.aX <- function(x, a, P, lambda, d, Ak, gamma.prev)
{
    if(!is.na(d))
    {
        # calculate these probabilities
        p.recomb.T <- P.recombination(TRUE, lambda, d)
        p.recomb.F <- 1 - p.recomb.T
    }else{
        p.recomb.T <- 1
        p.recomb.F <- 0
        gamma.prev[,] <- 0
    }

    # P(a | gamma)
    pag <- ifelse(is.na(a), 1, ifelse(a == 0, 1 - P[x], P[x]))

    # P(gamma)
    pg <- gamma.prev[,x] * p.recomb.F + Ak[,x] * p.recomb.T

    ## # P(a)
    ## pa <- matrix(0, nrow = length(a), ncol = length(P))

    ## # -- start with pag for each population
    ## for(l in 1:length(P))
    ## {
    ##     pa[,l] <- ifelse(a == 0, 1 - P[l], P[l])
    ## }

    ## # -- now multiply pg (Ak) through
    ## pa <- pa * Ak

    ## # sum up each row for pa
    ## pa <- apply(pa, 1, sum)

    return(pag * pg)# / pa)
}

P.a.gamma <- function(prev, geno, prior, P, Ak, lambda, d)
{
    pag <- array(NA, dim = c(dim(geno)[1], 2, length(P)),
                 dimnames = list(dimnames(geno)[[1]], c('Mother', 'Father'), names(P)))

    ### P(recombination with previous marker in MCMC chain) ###
    if(!is.na(d))
    {
        pr <- ## cbind(P.recombination(TRUE, lambda[,1], d),
                    P.recombination(TRUE, lambda[,2], d)#)
    }else{
        pr <- 1
    }

    pg <- pr * Ak + (1 - pr) * prev

    if(length(prior[[1]]$linked) > 0) # this shouldn't be possible, but somehow it is...???
    {
        pa <- pg # same format

        #### P(a) ####
        for(c in 1:2)
        {
            if(is.matrix(prior[[1]]$eig) | length(prior[[1]]$eig) > 1)
            {
                pcs <- geno[,prior[[1]]$linked,c] %*% prior[[1]]$eig
            }else{
                pcs <- geno[,prior[[1]]$linked,c]
            }

            lr <- cbind(1, pcs) %*% prior[[1]]$betas

            rat <- exp(lr)
            pa[,c,1] <- ifelse(rat == Inf, 1, exp(lr) / (1 + exp(lr)))
        }

        pa[,,2] <- 1 - pa[,,1]

        pga <- pa * pg

        for(c in 1:2)
            pga[,c,] <- pga[,c,] / apply(pga[,c,], 1, sum)

        return(pga)
    }else{
        return(pg)
    }
}

P.gammas.a <- function(locus, geno, prior, P, A0, Ak, lambda, d, gammas.prev)
{
    combos <- numeric()
    for(k1 in 1:length(P))
    {
        for(k2 in k1:length(P))
        {
            combos <- c(combos, paste('g', k1, k2, sep = ''))
        }
    }

    pga <- array(0, dim = c(dim(geno)[1], length(combos)),
                 dimnames = list(dimnames(geno)[[1]], combos))

    ### P(recombination with previous marker in MCMC chain) ###
    if(!is.na(d))
    {
        pr <- cbind(P.recombination(TRUE, lambda[,1], d),
                    P.recombination(TRUE, lambda[,2], d))

        # probability of one recombination
        tmp <- pr[,1] * (1 - pr[,2]) + (1 - pr[,1]) * pr[,2]

        # probability of two recombinations
        pr[,2] <- apply(pr, 1, prod)

        pr[,1] <- tmp
    }else{
        pr <- cbind(rep(0, dim(lambda)[1]), rep(1, dim(lambda)[1]))
    }

    ### P(O | known ancestral state) ###
    for(k in 1:length(combos))
    {
        # we will want these later
        k1 <- as.numeric(substr(combos[k], 2, 2))
        k2 <- as.numeric(substr(combos[k], 3, 3))

        # add probability of no recombinations on either chromosome
        pga[,k] <- (1 - apply(pr, 1, sum)) * gammas.prev[,k]

        # add probability of one recombination on one chromosome
        for(kprev in combos)
        {
            case <- c(combos[k] == kprev, # same
                      sum(strsplit(combos[k], '')[[1]] %in% strsplit(kprev, '')[[1]]) == 1) # couble recomb

            pga[,k] <- pga[,k] +
                switch(c('same', 'double', 'one')[c(case[1], case[2], all(!case))],
                       same = pr[,1] * ((A0[,k1] + A0[,k2]) / 2) * gammas.prev[,kprev],
                       double = 0,
                       one =
                       { # otherwise it was a specific one...which one did it change to?
                           changed <- c(k1, k2)[which(!c(k1, k2) %in% strsplit(kprev, '')[[1]] |
                                                      !strsplit(kprev, '')[[1]][-1] %in% c(k1, k2))]
                           pr[,1] * A0[,changed]  * gammas.prev[,kprev]
                       })
        }

        # add probability of double recombination to k
        pga[,k] <- pga[,k] + pr[,2] * (Ak[,1,k1]*Ak[,2,k2] + Ak[,1,k2]*Ak[,2,k1]) /
                                          (1 + (k1 == k2)) # divide by 2 when k1 == k2

        # probability of observing varant allele given ancestral state is k
        pO <- P.O2(geno[,prior$model[[k]]$linked, drop = FALSE], P = P[c(k1, k2)],
                   prior = prior$model[[k]], linked = length(prior$model[[k]]$linked) > 0)

        # " but with recombination it is no longer linked to previous markers
        pO2 <- P.O2(geno[,prior$model[[k]]$linked, drop = FALSE], P = P[c(k1, k2)],
                    prior = prior$model[[k]], linked = FALSE)

        # probability of a crossover within the block supporting pO
        p.b.recomb <- 1 - apply(exp(-lambda*prior$model[[k]]$d / 100), 1, prod)

        # probability of observed allele given ancestral state is k
        tmp <- cbind(locus == 0, locus == 1, locus == 2)

        pag <- ifelse(is.na(locus), 1, pO[tmp] * (1 - p.b.recomb) +
                                       pO2[tmp] * p.b.recomb)
        pga[,k] <- pga[,k] * pag
    }

    # normalize
    pga <- pga / apply(pga, 1, sum)

    return(pga)
}


###########################################
# P(n.crossovers == 0 ; lambda, distance) #
###########################################

P.0crossovers <- function(lambda, distance)
{
   exp(-lambda*distance / 100)
}

qodd.crossovers <- function(p.crossovers, lambda, distances)
{
    # we'll use these below
    p.recomb.T <- P.recombination(TRUE, lambda, distances)
    ld <- lambda*distances
    lld <- log(lambda*distances)

    retval <- rep(0, length(p.crossovers))

    for(j in 1:length(p.crossovers))
    {
        if(p.crossovers[j] >= 1)
        {
            retval[j] <- NA
            next
        }
        if(distances[j] == 0)
        {
            retval[j] <- 0
            next
        }

        p <- ld[j] * exp(-ld[j]) / p.recomb.T[j]

        n <- 3 # for next step -- always 1 step ahead here for the while loop
        tmp <- exp(n*lld[j] - ld[j] - lgamma(n)) / p.recomb.T[j]

        while(p + tmp <= p.crossovers[j])
        {
            p <- p + tmp

            n <- n + 2
            tmp <- exp(n*lld[j] - ld[j] - lgamma(n)) / p.recomb.T[j]
        }

        retval[j] <- n - 2 # we'll always go one step too far
    }

    return(retval)
}

##################################################################
# P(n.crossovers > 0 | gamma_j, gamma_j-1 ; Ak, lambda, distance #
##################################################################

# G = sampled ancestry states
# A = parent's global ancestry
# lambda = parent's lambda value
# distance = d -- should have 0's instead of NAs for 1st marker in each chromosome
# m = dim(G)[1] == length(distance)
P.1plus.crossover.G <- function(G, A, lambda, distance, m)
{
    # get population numbers for G
    g <- G %*% 1:dim(G)[2]

    # [P]robability of [0] and [P]robability of [1] [P]lus crossovers
    p0 <- exp(-lambda * distance[-1])
    p1p <- 1 - p0

    retval <- ifelse(g[-1] != g[-m], 1,
                     A[g[-1]] * p1p / (p0 + A[g[-1]] * p1p))

    return(retval)
}

P.1plus.crossover.gamma.gamma.prev <- function(gamma, gamma.prev, A, lambda, distance)
{
    # [P]robability of [0] and [P]robability of [1] [P]lus crossovers
    p0 <- P.0crossovers(lambda, distance)

    retval <- rep(0, dim(gamma)[1])

    for(g1 in 1:length(A))
        for(g2 in 1:length(A))
        {
            if(g1 != g2)
            {
                retval <- retval + gamma[,g1] * gamma.prev[,g2]
            }else{
                tmp <- exp(-2*lambda*distance)
                odd <- A[g1] * (1 - tmp) / 2
                even <- (1 + tmp) / 2 - p0
                retval <- retval + ( (odd + even) / (p0 + odd + even) ) * ( gamma[,g1] * gamma.prev[,g2] )
            }
        }

    return(retval)
}


############################
# P(lambdas | alpha, beta) #
############################

# lambdas = {lambda_11, lambda_12, lambda_21, lambda_22, ... , lambda_m1, lambda_m2}
logP.lambdas.hyper <- function(lambdas, alpha, beta, epsilon = 0.01)
{
    return(dgamma(lambdas, alpha, beta, log = TRUE))

    ## retval <- dgamma(lambdas, alpha, beta, log = TRUE)
    ## maxval <- pgamma(epsilon, alpha, beta, log = TRUE) # if we have 0's, this will avoid returning Inf

    ## return(ifelse(retval > maxval & lambdas < epsilon, maxval, retval))
}

##################
# P(P | tao, Pm) #
##################

# P in (0,1) = ancestral reference allele count
logP.P.tau.Pm <- function(P, tau, Pm)
{
    dbeta(P, tau*Pm, tau*(1 - Pm), log = TRUE)
}

##################
# P(Ajk | omega) #
##################

# Aks = matrix of local ancestry estimates
# Aks are dirichlet distributed with all parameters in alpha equal to omega
logP.Ak.omega <- function(Aks, omega, ratio = TRUE)
{
    # this cancels out in the likelihood ratio, since omega is constant
    if(ratio)
        return(log(Aks) %*% (omega - 1))

    return(lgamma(sum(omega)) - sum(lgamma(omega)) + log(Aks) %*% (omega - 1))
}

####################
# P(delta | sigma) #
####################

logP.delta.sigma <- function(deltas, sigma, lims, dev, epsilon = 0.01)
{
    if(length(sigma) == 1)
    {
        pdf <- ifelse(deltas > lims, 0, dnorm(deltas / sigma))
        return(log(pdf))
    }else{
        pdf <- ifelse(deltas > lims, 0, dnorm(deltas / rep(sigma, each = nrow(deltas))))

        if(is.null(dim(pdf)))
            return(log(pdf))

        return(apply(log(pdf), 1, sum))
    }
}

#####################################
# P(Ajk | omega; j on X chromosome) #
#####################################

# AX = matrix of local ancestry estimates
# Aks = matrix of local ancestry estimates
# Aks are dirichlet distributed with all parameters in alpha equal to omegaX*Aks
logP.AX.omega <- function(AX, Aks, omegaX, ratio = TRUE)
{
    if(ratio) # this cancels out in the likelihood ratio, since omegaX and Aks are constant
        return(apply(log(AX) * (omegaX*Aks - 1), 1, sum))

    return(lgamma(apply(omegaX*Aks, 1, sum)) - apply(apply(omegaX*Aks, 2, lgamma), 1, sum) +
           apply(log(AX) * (omegaX*Aks - 1), 1, sum))
}

###########################
# logP(gammas & geno | P) #
###########################

logP.gammas.geno.P <- function(gammas, geno, P)
{
    return(sapply(1:dim(P)[1], logP.gammas.geno.P.j, gammas = gammas, geno = geno, P = P))
}

logP.gammas.geno.P.j <- function(j, gammas, geno, P)
{
    # double check the type on gammas, too??? ...should always be numeric, right?
    .Call("logP_gammas_geno_P_j", gammas[,j,1,], gammas[,j,2,], geno[[j]], as.numeric(P[j,]))
}

####################################
# logP(gammas | lambda, distances) #
####################################

# gammas for one parent of one individual (one row per SNP, one column per population)
# single lambda value
# distance for all SNPs
logP.lambda.gammas <- function(gammas, lambda, distances)
{
    n <- length(distances)

    p.same <- apply(gammas[-1,] * gammas[-n,], 1, sum)
    p.diff <- 1 - p.same

    recombT <- P.recombination(TRUE, lambda, distances[-1])
    recombF <- P.recombination(FALSE, lambda, distances[-1])

    ps <- p.same*recombT + p.diff*recombF
    ps <- ps[distances[-1] > 0] # these are uninfrmative

    return(sum(log(ps), na.rm = TRUE))
}

########################
# rdirichlet(n, alpha) #
########################

# n - number of random variables to generate (i.e. number of rows)
# alpha - matrix of hyper-parameters governing dirichlet distribution (one row per random variable)

# !!! n should equal number of rows in alpha !!!
rdirichlet <- function(n, alpha)
{
    # random gammas
    tmp <- matrix(rgamma(length(alpha), as.vector(alpha)), nrow = n)

    # -> random dirichlet
    return(tmp / rowSums(tmp))
}


##########################################################################
# probability of allele state (variant) given gamma and previous markers #
##########################################################################

# This isn't for use in trying to assign alleles, but rather for identifying the
# likelihood of observing a variant allele under no recombination.

# This is current-marker-genotype-agnostic, but does rely on previous markers.

# If the marker is unlinked to previous markers, this simplifies to the values in P,
# otherwise it relies on linkage in the prior sample

P.O <- function(O.prev, P, prior, linked)
{
    dims <- dim(O.prev)
    if(length(dims) == 2)
        dims <- c(dims[1], 1, dims[2]) # if only one marker (or none), it will be mising it's second dimension

    .Call("P_O", as.numeric(O.prev), as.numeric(P), as.numeric(prior$betas),
          as.numeric(prior$eig), as.integer(linked & !is.null(prior$eig)),
          as.integer(dims))
}

P.O2 <- function(O.prev, P, prior, linked)
{
    pO <- matrix(NA, dim(O.prev)[1], 3, dimnames = list(dimnames(O.prev)[[1]], c('a=0', 'a=1', 'a=2')))

    # do this by individual... transfer to C
    if(linked & !is.null(prior$eig)) # sometimes we get uninformative markres nearby -> in pcr.prior
    {
        tmp <- cbind(O.prev[,] == 1, O.prev[,] == 2)

        pcs <- tmp %*% prior$eig
        lratios <- cbind(1, pcs) %*% prior$betas
        ratios <- exp(lratios)

        pO[,1] <- 1 / (1 + apply(ratios, 1, sum))
        pO[,2:3] <- ratios * pO[,1]

        if(!all(is.finite(pO)))
        {
            bad <- pO[,1] == 0
            pO[bad,'a=1'] <- 1 / (1 + exp(lratios[bad,1] - lratios[bad,2]))
            pO[bad,'a=2'] <- 1 / (1 + exp(lratios[bad,2] - lratios[bad,1]))
        }

        pO <- ifelse(pO == Inf, 1, pO)
    }else{
        pO[,'a=0'] <- (1 - P[1]) * (1 - P[2])
        pO[,'a=1'] <- P[1] * (1 - P[2]) + (1 - P[1]) * P[2]
        pO[,'a=2'] <- P[1] * P[2]
    }

    return(pO)
}

################################
# Inverse Wishart distribution #
################################

# Based on code from MCMCpack
rwish <- function(v, S, Sinv = NULL)
{
    # not many checks...this assumes it has been given correctly formatted data objects
    p <- nrow(S)
    if(v < p)
        v <- p

    # this is actually Wishart, so don't forget to invert for our needs
    if(is.null(Sinv))
        Sinv <- solve.approx(S)

    CC <- try(chol(Sinv), silent = TRUE)

    if(class(CC) == 'try-error')
        return(rwish2(v, Sinv))

    Z <- matrix(0, p, p)

    diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
    if(p > 1) {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }

    return(crossprod(Z %*% CC))
}

# adapted from code in a blog post by Simon Barthelme
rwish2 <- function(v,Sinv)
{
    X <- rmvnorm(v,sig=Sinv)
    return(t(X)%*%X)
}

# Based on code from MCMCpack
# log scale
ldiwish <- function(W, Winv, v, S)
{
    # not many checks...this assumes it has been given correctly formatted data objects
    k <- nrow(S)
    if(v < k)
        v <- k

    # denominator
    gammapart <- 0
    for(i in 1:k) {
      gammapart <- gammapart + lgamma((v + 1 - i)/2)
    }
    denom <- gammapart + (v * k / 2) * log(2) + (k*(k-1)/4) * log(pi)

    # numerator
    detS <- det(S)
    detW <- det(W)
    hold <- S %*% Winv
    tracehold <- sum(hold[row(hold) == col(hold)])
    num <- (v/2) * log(detS) + (-(v + k + 1)/2)*log(detW) + (-1/2 * tracehold)

    return(num - denom)
}


#########################################################
# Likelihood ratio comparing new betas with prior betas #
#########################################################

logP.beta.sigma.tau <- function(beta, beta.p, sigma, tau)
{
    return(.Call("logP_beta_sigma_tauR", beta, beta.p, sigma, tau))

    # (2pi)^-(length(beta)/2) * det(Sigma)^(-1/2)
    # will cancel out in the ratio, so we don't include that here...
    ## return((length(beta) * log(tau) - sum(log(sigma) + tau * (beta - beta.p)^2 / sigma)) / 2)
}


#####################################################################
# Log Likelihood of sampled ancestral states, G, given lambda and d #
#####################################################################

logP.G.lambda.d <- function(G, Ak, lambda, d, phased)
{
    retval <- G[,1,,1]
    retval[] <- 0


    for(i in 1:dim(G)[1])
    {
        lambda.i <- mean(lambda[i,], na.rm = TRUE)
        A.i <- apply(Ak[i,,], 2, mean)

        for(k in 1:dim(G)[3])
        {
            loglik <- 0

            for(j in 2:dim(G)[2])
            {
                if(phased) # when we have phased data
                {
                    P.r <- P.recombination(TRUE, lambda[i,k], d[j])

                    loglik <- loglik + log(P.r * Ak[i,k,which(as.logical(G[i,j,k,]))] + (1 - P.r))

                    if(any(G[i,j-1,k,] != G[i,j,k,]))
                    {
                        retval[i,k] <- retval[i,k] + loglik + log(1 - exp(loglik))

                        loglik <- 0
                    }
                }else{ # when we have genotype data (unphased)
                    P.r <- P.recombination(TRUE, lambda.i, d[j])

                    anc <- as.numeric(strsplit(dimnames(G)[[4]][which(as.logical(G[i,j,1,]))], '')[[1]][2:3])

                    loglik <- loglik + log(P.r * A.i[anc[1]] + (1 - P.r)) + log(P.r * A.i[anc[2]] + (1 - P.r))

                    if(any(G[i,j-1,k,] != G[i,j,k,]))
                    {
                        retval[i] <- retval[i] + loglik + log(1 - exp(loglik))

                        loglik <- 0
                    }
                }

            }
        }
    }

    return(retval)
}
