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


###############################################################################################
# P(gamma == x | a ; P, lambda, lambda.alt, distance, Ak, Ak.alt, gamma.prev, gamma.alt.prev) #
###############################################################################################

P.gammas <- function(geno, P, Pm.prior, lambda, lambdaX, d, chr, A0, Ak, AX, gender, sex.chr, unphased, dev)
{
    if(!unphased)
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
            ### P(gammas | a) ###
            gammas[,j,,] <- P.gamma.a(gammas[,max(j-1, 1),,], geno, Pm.prior[[j]],
                                      P[j,], Ak, lambda, d[j], unphased)
        }

        g.prev <- gammas[,1,,] # ignored, but this is the correct format...

        # Reverse Chain #
        for(j in ord[length(ord):1])
        {
            g.prev <- P.a.gamma(g.prev, geno, Pm.prior[[j]], P[j,], Ak, lambda, d[j+1], unphased)
            gammas[,j,,] <- gammas[,j,,] * g.prev
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

P.gamma.a <- function(prev, geno, prior, P, Ak, lambda, d, unphased)
{
    ### P(recombination with previous marker in MCMC chain) ###
    if(!is.na(d))
    {
        pr <- ## cbind(P.recombination(TRUE, lambda[,1], d),
                    P.recombination(TRUE, lambda[,2], d)#)
    }else{
        pr <- 1
    }

    if(unphased)
    {
        # r = 0
        pg <- (1 - pr)^2 * prev

        for(k1 in 1:dim(Ak)[3])
        {
            for(k2 in k1:dim(Ak)[3])
            {
                index <- paste('g', k1, k2, sep = '')

                # r = 1
                p1 <- prev[,substr(dimnames(prev)[[2]], 2,2) == k1, drop = FALSE]
                p2 <- prev[,substr(dimnames(prev)[[2]], 3,3) == k2, drop = FALSE]

                if(dim(p1)[2] > 0)
                    pg[,index] <- pg[,index] + 2*pr*(1 - pr) * apply(p1 * Ak[,1,k1] / 2, 1, sum)

                if(dim(p2)[2] > 0)
                    pg[,index] <- pg[,index] + 2*pr*(1 - pr) * apply(p2 * Ak[,1,k2] / 2, 1, sum)

                # r = 2
                pg[,index] <- pg[,index] + pr^2 * Ak[,1,k1] * Ak[,1,k2] * ifelse(k1 == k2, 1, 2) # r = 2
            }
        }
    }else{
        pg <- pr * Ak + (1 - pr) * prev
    }

    if(length(prior$linked) > 0) # this shouldn't be possible, but somehow it is...???
    {
        pag <- pg # same format

        #### P(a) ####
        index <- array(FALSE, dim = dim(geno), dimnames = dimnames(geno))
        if(unphased)
        {
            index[,prior$linked] <- TRUE
        }else{
            index[,prior$linked,] <- TRUE
        }

        # calculate PCR predictors based on PC loadings
        if(is.matrix(prior$eig) | length(prior$eig) > 1)
        {
            pcs <- geno[index] %*% prior$eig
        }else{
            pcs <- geno[index]
        }

        # log likelihood ratio
        lr <- cbind(1, pcs) %*% prior$betas

        # ratio
        rat <- exp(lr)

        # probability of a given g
        index <- array(TRUE, dim = dim(pag), dimnames = dimnames(pag))
        if(unphased)
        {
            index[,dim(index)[2]] <- FALSE
        }else{
            index[,,dim(index)[3]] <- FALSE
        }

        denom <- apply(lr, 1, sum)
        pag[index] <- ifelse(rat == Inf, 1, exp(lr) / (1 + exp(denom)))
        pag[!index] <- 0

        pag[!index] <- 1 - apply(pag, 1:(length(dim(pag)) - 1), sum)

        pga <- pag * pg

        # divide by sum over the last dimension
        pga <- pga / apply(pga, 1:(length(dim(pga)) - 1), sum)

        return(pga)
    }else{
        return(pg)
    }
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
