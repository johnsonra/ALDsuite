# updates.R
# Gibbs samplers to update parameters during the Monte-Carlo step
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created September 5, 2012
# Last Modified March 24, 2014

#########################################################################
# Gibbs sampler for updating alpha and beta - defining prior for lambda #
#########################################################################

# alpha = study-wide estimate of alpha
# beta = study-wide estimate of beta
# lambdas = vector of values for lambda for all parents (one row per individual, one column per parent)
# step = parameter for scaling the step size between guesses
# burn = number of burn-in iterations
# iter = number of follow-on iterations
# every = save results afer each "every" cycles
# debug = indicator used to return extra information when debugging this function
updt.alpha.beta <- function(alpha, beta, lambdas, step = 1, burn = 50, iter = 100,
                            every = 10, debug = FALSE)
{
  # optimizing in terms of mean and variance
  m <- alpha / beta
  v <- alpha / beta^2

  # current log-likelihood
  loglik.cur <- sum(logP.lambdas.hyper(lambdas, m^2/v, m/v))

  # burn in
  burn.m <- numeric(floor(burn / every))
  burn.v <- numeric(floor(burn / every))

  # follow on
  iter.m <- numeric(floor(iter / every))
  iter.v <- numeric(floor(iter / every))

  # run all cycles
  for(i in 1:(burn + iter))
  {
    # update m
    m.new <- 10^rnorm(1, mean = 1, sd = 0.5)
    loglik.new <- sum(logP.lambdas.hyper(lambdas, m.new^2/v, m.new/v))

    p.switch <- min(c(1, exp(loglik.new - loglik.cur)))

    if(is.na(p.switch))
    {
        save(list = ls(envir = environment()), file = 'bad_m.new.RData', envir = environment())
        stop('Found a bad m.new in updt.alpha.beta()')
    }

    if(!is.na(p.switch) & p.switch >= 0 & p.switch <= 1)
        if(rbinom(1, 1, p.switch))
        {
            m <- m.new
            loglik.cur <- loglik.new
        }

    # update v
    v.new <- 10^rnorm(1, mean = 1, sd = 0.5)
    loglik.new <- logP.lambdas.hyper(lambdas, m^2/v.new, m/v.new)

    p.switch <- min(c(1, exp(loglik.new - loglik.cur)))

    if(is.na(p.switch))
    {
        save.image(file = 'bad_v.new.RData')
        stop('Found a bad v.new in updt.alpha.beta()')
    }

    if(!is.na(p.switch) & p.switch >= 0 & p.switch <= 1)
        if(rbinom(1, 1, p.switch))
        {
            v <- v.new
            loglik.cur <- loglik.new
        }

    # save on every
    if(i <= burn)
    {
      if(i / every == floor(i / every))
      {
        burn.m[i / every] <- m
        burn.v[i / every] <- v
      }
    }else{
      if((i - burn) / every == floor((i - burn) / every)){
        iter.m[(i - burn) / every] <- m
        iter.v[(i - burn) / every] <- v
      }
    }
  }

  if(debug)
    return(list(burn.m = burn.m, burn.v = burn.v, iter.m = iter.m, iter.v = iter.v))
  else
    return(list(alpha = mean(iter.m)^2 / mean(iter.v), beta = mean(iter.m) / mean(iter.v)))
}


#######################################################
# Update of individual lambdas (one lambda at a time) #
#######################################################

# probably want to remove some variables here (G0, A0, sigma...P.Aeq1?)
updt.lambda <- function(G, Ak, lambda, Pm.prior, distances, alpha, beta, chr, phased, dev, W = NULL) #### X chromosome???
{
    # if W is null, we are doing the autosomes
    if(is.null(W)){
        W <- rep(TRUE, dim(lambda)[1])
    }else{
        W <- 1:dim(lambda)[1] %in% W # these are the females
    }

    # include all individuals for mothers' lambdas
    include <- rep(TRUE, dim(lambda)[1])

    for(k in 1:2)
    {
        if(k == 2 & !phased)
        {
            lambda[include, 2] <- lambda[include, 1]
            next
        }

        # only include those specified by W for fathers' lambdas
        # (males don't have an X from fathers)
        if(k == 2)
            include <- W

        if(all(!W)) # just in case we have all male sex chromosomes
            next

        loglik.cur <- logP.lambdas.hyper(lambda[include, k], alpha, beta)

        lambda.new <- sample.lambda(k, G[include,,,, drop = FALSE], Ak[include,,], lambda[include,],
                                    Pm.prior, distances, chr, alpha, beta, phased, dev)

        # make decision to keep new or old estimate
        loglik.new <- logP.lambdas.hyper(lambda.new, alpha, beta)

        # make decision for each individual
        p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
        b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

        # update
        lambda[include,k][b.switch] <- lambda.new[b.switch]
    }

    return(lambda)
}

################################################################################################
# Update sample of P with respect to the modern frequencies, Pm, and the dispersion factor Tau #
################################################################################################

## # P = matrix of ancestral reference allele frequencies (one row per SNP, one column per population)
## # tau = vector of dispersion factors (one for each ancestral population)
## # Os = matrix of observed reference alleles (one row for each individual, one column for each population)
## # cl = cluster object for parallelization
## updt.P <- function(P, tau, Pm, Os, coolness, dev)
## {
##     P.new <- P

##     # sum up "observed" alleles
##     ref <- apply(Os[,,,1], 2:3, sum)
##     vart <- apply(Os[,,,2], 2:3, sum)

##     # sample new value of P
##     for(l in 1:dim(P)[2])
##     {
##         alpha <- Pm[,l]*tau[l] + ref[,l] + 0.5
##         beta <- (1 - Pm[,l])*tau[l] + vart[,l] + 0.5

##         P.new[,l] <- rbeta(dim(P)[1], alpha, beta)

##         loglik.cur <- logP.P.tau.Pm(P[,l], tau[l], Pm[,l])
##         loglik.new <- logP.P.tau.Pm(P.new[,l], tau[l], Pm[,l])

##         # decide
##         p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
##         b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch)) # decide on each one

##         # update -- weight by coolness
##         P[b.switch,] <- P[b.switch,] * (1 - coolness) + P.new[b.switch,] * coolness
##     }

##     return(P)
## }

updt.P.betas <- function(P, Pm, Pm.prior, lambda, tau, G, haps, phased, dev)
{
    # sample new value of P that are unlinked in either direction
    # there will be atleast 1 unlinked marker at each end of each chromosome
    for(k in 1:dim(P)[2])
    {
        P.old <- P[,k]
        Pm.k <- Pm[,k]

        # continuity correction
        ref <- rep(0.5, length(Pm.prior))
        vart <- ref

        # sum up observed alleles
        for(j in 1:length(Pm.prior))
        {
            if(phased)
            {
                ref[j] <- ref[j] + sum(1 - haps[,j,][as.logical(G[,j,,k])])
                vart[j] <- vart[j] + sum(haps[,j,][as.logical(G[,j,,k])])
            }else{
                for(combo in grep(k, dimnames(G)[[4]], value = TRUE))
                {
                    k1 <- as.numeric(substr(combo, 2, 2))
                    k2 <- as.numeric(substr(combo, 3, 3))

                    if(k1 == k2) # homozygous for k
                    {
                        ref[j] <- ref[j] + sum(2 - haps[as.logical(G[,j,1,combo]),j])
                        vart[j] <- vart[j] + sum(haps[as.logical(G[,j,1,combo]),j])
                    }else{ # heterozygous for k
                        hets <- sum(haps[as.logical(G[,j,1,combo]),j] == 1)

                        ref[j] <- ref[j] + sum(haps[as.logical(G[,j,1,combo]),j] == 0) +
                                           hets * (1 - P[j,k] / (P[j,k1] + P[j,k2]))
                        vart[j] <- vart[j] + sum(haps[as.logical(G[,j,1,combo]),j] == 2) +
                                             hets * P[j,k] / (P[j,k1] + P[j,k2])
                    }
                }
            }
        }

        alpha <- Pm.k*tau[k] + vart
        beta <- (1 - Pm.k)*tau[k] + ref

        P.new <- rbeta(length(P.old), alpha, beta)

        loglik.cur <- logP.P.tau.Pm(P.old, tau[k], Pm.k)
        loglik.new <- logP.P.tau.Pm(P.new, tau[k], Pm.k)

        # decide
        p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
        b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch)) # decide on each one

        # update -- weight by coolness
        P[b.switch,k] <- P.new[b.switch]
    }

    #set.seed(298347)
    # sample new values for betas
    for(j in 1:length(Pm.prior))
    {
        for(k in 1:length(Pm.prior[[j]]$model))
        {
            if(phased)
            {
                tmp <- tau[k]
            }else{
                k1 <- as.numeric(substr(dimnames(G)[[4]][k], 2, 2))
                k2 <- as.numeric(substr(dimnames(G)[[4]][k], 3, 3))
                tmp <- mean(tau[c(k1, k2)])
            }

            if(length(Pm.prior[[j]]$model[[k]]$linked) > 0)
                Pm.prior[[j]]$model[[k]] <- updt.betas(j, names(Pm.prior)[j], Pm.prior[[j]]$model[[k]],
                                                       lambda, G, haps, k, tmp, phased, dev = dev)
        }
    }

    return(list(P = P,
                Pm.prior = Pm.prior))
}

updt.betas <- function(j, snp, prior, lambda, G, haps, k, tau, phased, dev)
{
    # anchor SNPs coverd by haplotype
    covered <- c(j, which(dimnames(G)[[2]] %in% prior$linked))

    # collect haps that are entirely from population k
    if(phased)
    {
        moms2keep <- as.logical(apply(G[, covered, 'Mother', k, drop = FALSE], 1, prod))
        dads2keep <- as.logical(apply(G[, covered, 'Father', k, drop = FALSE], 1, prod))

        n <- sum(moms2keep) + sum(dads2keep)

        obs.haps <- rbind(haps[moms2keep, c(prior$linked, snp), 1],
                          haps[dads2keep, c(prior$linked, snp), 2])

        # max of 1 and (lower bound of window - 1)
        tmp <- max(min(covered) - 1, 1)
        need.wts <- !as.logical(c(G[, tmp, 'Mother', k] * moms2keep, G[, tmp, 'Father', k] * dads2keep))

        # min of number of anchor markers and (upper bound of window + 1)
        tmp <- min(max(covered) + 1, dim(G)[2])
        needs.wts <- !as.logical(c(G[, tmp, 'Mother', k] * moms2keep, G[, tmp, 'Father', k] * dads2keep))

        # calculate weights based on probability of recombination within the window covering the haplotype
        # assign 1 for haps where G is the same at markers neighboring the window
        # assign prob of crossover otherwise
        wts <- ifelse(needs.wts[c(moms2keep, dads2keep)],
                      ppois(0, prior$d * c(lambda[,1], lambda[,2]) / 100), 1)

        # get estimates for current sample
        if(length(wts) > 0)
        {
            if(length(prior$linked) == 1)
            {
                model <- logitreg(obs.haps[,snp], matrix(obs.haps[,prior$linked]), wt = wts, hessian = FALSE)
            }else{
                pcs <- obs.haps[,prior$linked] %*% prior$eig
                model <- logitreg(obs.haps[,snp], pcs, wts, hessian = FALSE)
            }
        }else{
            model <- list(par = 0)
        }

    }else{
        tokeep <- as.logical(apply(G[,covered,,k, drop = FALSE], 1, prod))
        obs.haps <- haps[tokeep, c(prior$linked, snp), drop = FALSE]

        n <- sum(tokeep)

        # get estimates for current sample
        if(n > 1)
        {
            if(length(prior$linked) == 1)
            {
                model <- multilogitreg(cbind(obs.haps[,snp] == 0, obs.haps[,snp] == 1, obs.haps[,snp] == 2),
                                       cbind(obs.haps[,prior$linked] == 1, obs.haps[,prior$linked] == 2),
                                       hessian = FALSE)
            }else{
                pcs <- cbind(obs.haps[,prior$linked] == 1, obs.haps[,prior$linked] == 2) %*% prior$eig
                model <- multilogitreg(cbind(obs.haps[,snp] == 0, obs.haps[,snp] == 1, obs.haps[,snp] == 2),
                                       pcs, hessian = FALSE)
            }
        }else{
            model <- list(par = 0)
            n <- 0
        }
    }

    # sample value for beta
    betas.new <- rnorm(length(prior$betas), (n*model$par + tau*prior$betas.p) / (n + tau),
                       prior$sigma / (n + tau)^2)

    if(!phased)
    {
        betas.new <- matrix(betas.new, ncol = 2)
    }

    # make decision
    loglik.cur <- logP.beta.sigma.tau(prior$betas, prior$betas.p, prior$sigma, tau)
    loglik.new <- logP.beta.sigma.tau(betas.new, prior$betas.p, prior$sigma, tau)

    p.switch <- min(1, exp(loglik.new - loglik.cur))

    if(as.logical(rbinom(1, 1, p.switch)))
       prior$betas <- betas.new

    return(prior)
}

##############################################
# Gibbs sample of Tau and Pm  - priors for P #
##############################################

updt.tau.Pm.prep <- function(tau, Pm.prior)
{
    # these are the betas we're going to update
    to.update <- list()

    for(k in 1:length(Pm.prior[[1]]$model))
    {
        to.update[[k]] <- list(j = numeric(),
                               priors = list())

        for(j in 1:length(Pm.prior))
        {
            if(length(Pm.prior[[j]]$model[[k]]$linked) > 0 & !is.null(Pm.prior[[j]]$model[[k]]$eig))
            {
                to.update[[k]]$j <- c(to.update[[k]]$j, j)
                jj <- length(to.update[[k]]$j)

                to.update[[k]]$priors[[jj]] <- list()
                to.update[[k]]$priors[[jj]]$betas.p <- rep(0, length(Pm.prior[[j]]$model[[k]]$betas.p))
            }
        }

        names(to.update[[k]]$priors) <- names(Pm.prior)[to.update[[k]]$j]
    }

    return(to.update)
}

# tau = population dispersion factor for one population
# Pm = vector of modern reference allele frequencies for each population (doing one population at a time here!!)
# Pm.counts = matrix of modern allele counts (first column is reference, second column is variant)
# P = vector of ancestral reference allele frequencies
# burn = number of burn in cycles
# iter = number of follow on cycles
# every = save results afer each "every" cycles
# debug = indicator to return extra information when debugging this function
updt.tau.Pm <- function(tau, Pm, Pm.counts, Pm.prior, P, to.update, burn = 50, iter = 100,
                        every = 10, debug = FALSE, dev = FALSE)
{
    # 1's and 0's will cause problems at boundaries
    P <- ifelse(P == 0, 1e-10, P)
    P <- ifelse(P == 1, 1 - 1e-10, P)

    if(debug)
    {
        # burn in
        burn.Pm <- array(0, dim = c(dim(Pm), floor(burn / every)))
        burn.tau <- matrix(NA, nrow = length(tau), ncol = floor(burn / every))

        # follow on
        iter.Pm <- array(0, dim = c(dim(Pm), floor(iter / every)))
        iter.tau <- matrix(NA, nrow = length(tau), ncol = floor(iter / every))
    }

    # these are the updated parameters we will return
    Pm.out <- matrix(0, ncol = ncol(Pm), nrow = nrow(Pm), dimnames = dimnames(Pm))
    tau.out <- rep(0, length(tau))

    # run all cycles
    for(i in 1:(burn + iter))
    {
        for(k in 1:length(tau))
        {
            # current log-likelihood
            loglik.cur <- logP.P.tau.Pm(P[,k], tau[k], Pm[,k])

            ### update Pm ###
            Pm.new <- rbeta(dim(Pm.counts)[1], Pm.counts[,k,'vart'],
                                               Pm.counts[,k,'ref'])

            # the likelihood for each one will depend only on that observation -- all others will cancel out
            # -> may as well do them all at once.
            loglik.new <- logP.P.tau.Pm(P[,k], tau[k], Pm.new)

            p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
            b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch)) # decide on each one

            Pm[b.switch,k] <- Pm.new[b.switch]
            loglik.cur[b.switch] <- loglik.new[b.switch]

            # sample the new value for tau so we can collect LRs as we go along through the update of beta
            loglik.tau.cur <- sum(loglik.cur)

            tau.new <- 10^rnorm(1, mean = 2, sd = 0.5)

            loglik.tau.new <- sum(logP.P.tau.Pm(P[,k], tau.new, Pm[,k]))

            point <- 0
            ### update beta.p ###
            ## if(length(to.update[[k]]$j) > 0)
            ## {
            ##     for(j in to.update[[k]]$j)
            ##     {
            ##         jj <- which(to.update[[k]]$j == j)

            ##         tmp <- sample.tau.Pm.prior(Pm.prior[[j]]$model[[k]], tau[k], tau.new, dev)

            ##         Pm.prior[[j]]$model[[k]] <- tmp$prior

            ##         loglik.tau.cur <- loglik.tau.cur + tmp$loglik.tau.cur
            ##         loglik.tau.new <- loglik.tau.new + tmp$loglik.tau.new
            ##     }
            ## }

            ### update tau ###
            p.switch <- min(c(1, exp(loglik.tau.new - loglik.tau.cur)))

            if(rbinom(1, 1, p.switch))
                tau[k] <- tau.new
        }

        # if in debug mode, save things as we go along
        if(debug)
        {
            if(i <= burn)
            {
                if(i / every == floor(i / every))
                {
                    burn.Pm[,,i / every] <- Pm
                    burn.tau[,i / every] <- tau
                }
            }else{
                if((i - burn) / every == floor((i - burn) / every)){
                    iter.Pm[,,(i - burn) / every] <- Pm
                    iter.tau[,(i - burn) / every] <- tau
                }
            }
        }

        # save current sample on every
        if(i > burn & (i - burn) / every == floor((i - burn) / every))
        {
            Pm.out <- Pm.out + Pm / (iter / every)
            tau.out <- tau.out + tau / (iter / every)

            ## for(k in 1:length(tau))
            ## {
            ##     if(length(to.update[[k]]$j) > 0)
            ##     {
            ##         for(j in 1:length(to.update[[k]]$j))
            ##         {
            ##             jj <- to.update[[k]]$j[j]
            ##             to.update[[k]]$priors[[j]]$betas.p <- to.update[[k]]$prior[[j]]$betas.p +
            ##                                                   Pm.prior[[jj]]$model[[k]]$betas.p / (iter / every)
            ##         }
            ##     }
            ## }
        }
    }

    ## for(k in 1:length(tau))
    ##     if(length(to.update[[k]]$j) > 0)
    ##         for(j in 1:length(to.update[[k]]$j))
    ##             Pm.prior[[to.update[[k]]$j[j]]]$model[[k]]$betas.p <- to.update[[k]]$priors[[j]]$betas.p

    if(debug)
        return(list(burn.Pm = burn.Pm, burn.tau = burn.tau, iter.Pm = iter.Pm, iter.tau = iter.tau,
                    Pm = Pm.out, tau = tau.out))
    else
        return(list(Pm = Pm.out, tau = tau.out, Pm.prior = Pm.prior))
}

sample.tau.Pm.prior <- function(prior, tau, tau.new, dev)
{
    # segfault error encountered here at some point?
    ## retval <- .Call('sample_tau_Pm_prior', prior$betas, prior$betas.p, prior$betas.prior,
    ##                 prior$sigma, tau, tau.new)

    ## prior$betas.p <- retval[-(1:2)]
    ## return(list(prior = prior,
    ##             loglik.tau.cur = retval[1],
    ##             loglik.tau.new = retval[2]))

    loglik.cur <- logP.beta.sigma.tau(prior$betas, prior$betas.p, prior$sigma, tau)

    # new beta.p
    beta.new <- rnorm(length(prior$betas.p), prior$betas.prior, prior$sigma)

    loglik.new <- logP.beta.sigma.tau(prior$betas, beta.new, prior$sigma, tau)

    p.switch <- min(c(1, exp(loglik.new - loglik.cur)))

    if(as.logical(rbinom(1, 1, p.switch)))
    {
        prior$betas.p <- beta.new
        loglik.cur <- loglik.new
    }

    return(list(prior = prior,
                loglik.tau.cur = loglik.cur,
                loglik.tau.new = logP.beta.sigma.tau(prior$betas, prior$betas.p, prior$sigma, tau.new)))
}

#######################################################################
# Gibbs sampler for updating omega - hyper parameter for prior on Aks #
#######################################################################

updt.omega <- function(omega, A0, burn = 50, iter = 100, every = 10, debug = FALSE)
{
    # force A0 to remain within (0,1)
    A0s <- ifelse(A0 < 1e-10, 1e-10, A0)

    loglik.cur <- sum(logP.Ak.omega(A0s, omega, ratio = FALSE))

    # burn in
    burn.omega <- matrix(NA, floor(burn / every), length(omega))

    # follow on
    iter.omega <- matrix(NA, floor(iter / every), length(omega))

    # run all cycles
    for(i in 1:(burn + iter))
    {
        omega.new <- 10^rnorm(length(omega), mean = 1, sd = 0.5)
        omega.new <- ifelse(omega.new < 1, 1 + omega.new, omega.new) # don't go down this road!

        ## for(l in 1:length(omega))
        ## {
        ##     old <- rep(1, length(omega)); old[l] <- 0
        ##     new <- rep(0, length(omega)); new[l] <- 1

            ## loglik.new <- logP.Ak.omega(Aks, omega * old + omega.new * new, ratio = FALSE)
        loglik.new <- sum(logP.Ak.omega(A0s, omega.new, ratio = FALSE))

            p.switch <- min(c(1, exp(loglik.new - loglik.cur)))

            # make them equal and move to the next parameter
            if(rbinom(1, 1, p.switch)){
                ## omega[l] <- omega.new[l]
                omega <- omega.new
                loglik.cur <- loglik.new
            }
        ## }

        # save on every
        if(i <= burn)
        {
            if(i / every == floor(i / every))
                burn.omega[i / every,] <- omega
        }else{
            if((i - burn) / every == floor((i - burn) / every)){
                iter.omega[(i - burn) / every,] <- omega
            }
        }
    }

    if(debug)
        return(list(omega = apply(iter.omega, 2, mean), burn = burn.omega, iter = iter.omega))
    else
        return(apply(iter.omega, 2, mean))
}

updt.omegaX <- function(omegaX, Ak, AX, M, burn = 50, iter = 100, every = 10, debug = FALSE)
{
    # unlist Ajk
    Aks <- apply(Ak, c(1,3), mean)
    Aks[M,] <- Ak[M,1,]
    Aks <- ifelse(Aks < 1e-10, 1e-10, Aks)

    AXs <- apply(AX, c(1,3), mean, na.rm = TRUE)
    AXs <- ifelse(AXs < 1e-10, 1e-10, AXs)

    loglik.cur <- sum(logP.AX.omega(AXs, Aks, omegaX, ratio = FALSE))

    # burn in
    burn.omegaX <- numeric(floor(burn / every))

    # follow on
    iter.omegaX <- numeric(floor(iter / every))

    # run all cycles
    for(i in 1:(burn + iter))
    {
        omegaX.new <- 10^rnorm(1, mean = 1, sd = 0.5)

        loglik.new <- sum(logP.AX.omega(AXs, Aks, omegaX.new))

        p.switch <- min(c(1, exp(loglik.new - loglik.cur)))

        # make them equal and move to the next parameter
        if(rbinom(1, 1, p.switch)){
            omegaX <- omegaX.new
            loglik.cur <- loglik.new
        }


        # save on every
        if(i <= burn)
        {
            if(i / every == floor(i / every))
                burn.omegaX[i / every] <- omegaX
        }else{
            if((i - burn) / every == floor((i - burn) / every)){
                iter.omegaX[(i - burn) / every] <- omegaX
            }
        }
    }

    if(debug)
        return(list(omegaX = mean(iter.omegaX), burn = burn.omegaX, iter = iter.omegaX))
    else
        return(mean(iter.omegaX))
}


###################################################################################
# Update A0, given all local ancestry estimates and distances between the markers #
###################################################################################

updt.A <- function(A0, Ak, G, omega, sigma, epsilon, phased, dev)
{
    loglik.cur <- logP.Ak.omega(A0, omega)

    if(phased)
    {
        ### update A ###

        totals1 <- apply(G[,,1,], c(1,3), sum)  # summing over markers for each parent
        totals2 <- apply(G[,,2,], c(1,3), sum)

            # sample new candidate for A
        A1.new <- rdirichlet(dim(A0)[1], omega + totals1 + 0.5)
        A2.new <- rdirichlet(dim(A0)[1], omega + totals2 + 0.5)

        A.new <- (A1.new + A2.new) / 2

        # make decision to keep new or old estimate
        loglik.new <- logP.Ak.omega(A.new, omega)

        p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
        b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

        # finalize update of A0
        A0[b.switch,] <- A.new[b.switch,]

        Ak[b.switch,1,] <- A1.new[b.switch,]
        Ak[b.switch,2,] <- A2.new[b.switch,]

        ### probability of no admixture ###
        P.Aeq1 <- Ak # same dimensions
        P.Aeq1[,,] <- NA

        # this will be beta distributed for each individual population (e.g. African vs something else)
        if(dim(P.Aeq1)[3] == 2)
        {
            for(l in 1:2)
            {
                P.Aeq1[,1,l] <- pbeta(1 - epsilon, omega[l] + totals1[,l],
                                      omega[-l] + totals1[,-l], lower.tail = FALSE)
                P.Aeq1[,2,l] <- pbeta(1 - epsilon, omega[l] + totals2[,l],
                                      omega[-l] + totals2[,-l], lower.tail = FALSE)
            }
        }else{
            for(l in 1:dim(P.Aeq1)[3])
            {
                P.Aeq1[,1,l] <- pbeta(1 - epsilon, omega[l] + totals1[,l],
                                      sum(omega[-l]) + apply(totals1[,-l], 1, sum), lower.tail = FALSE)
                P.Aeq1[,2,l] <- pbeta(1 - epsilon, omega[l] + totals2[,l],
                                      sum(omega[-l]) + apply(totals2[,-l], 1, sum), lower.tail = FALSE)
            }
        }
    }else{
        ##### Update A for unphased data #####
        # don't worry about Ak for now...
        totals <- matrix(0, nrow = dim(G)[1], ncol = length(omega))

        for(k1 in 1:length(omega))
        {
            for(k2 in k1:length(omega))
            {
                tmp <- apply(G[,,paste('g', k1, k2, sep = '')], 1, sum)
                totals[,k1] <- totals[,k1] + tmp
                totals[,k2] <- totals[,k2] + tmp
            }
        }

        # new sample
        A.new <- rdirichlet(dim(A0)[1], omega + totals + 0.5)

        # updated likelihood
        loglik.new <- logP.Ak.omega(A.new, omega)

        # make decision
        p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
        b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

        A0[b.switch,] <- A.new[b.switch,]

        # updates of Ak...not tracking these separately for now...
        Ak[,1,] <- A0
        Ak[,2,] <- A0

        # P(A == 1)
        P.Aeq1 <- Ak # same dims
        P.Aeq1[,,] <- NA

        for(k in 1:dim(P.Aeq1)[3])
        {
            P.Aeq1[,1,k] <- pbeta(1 - epsilon, omega[k] + totals[,k],
                                  sum(omega[-k]) + apply(totals[,-k,drop = FALSE], 1, sum),
                                  lower.tail = FALSE)
        }

        P.Aeq1[,2,] <- P.Aeq1[,1,]
    }

    return(list(A0 = A0, Ak = Ak, sigma = sigma, P.Aeq1 = P.Aeq1))
}

## updt.A <- function(A0, Ak, Ajs, omega, sigma, epsilon, dev)
## {
##     loglik.cur <- logP.Ak.omega(A0, omega)

##     ### update A ###

##     totals1 <- apply(Ajs[,,1,], c(1,3), sum)  # summing over markers for each parent
##     totals2 <- apply(Ajs[,,2,], c(1,3), sum)
##     totals <- totals1 + totals2

##     # sample new candidate for A
##     A.new <- rdirichlet(dim(A0)[1], omega + totals + 0.5)

##     # make decision to keep new or old estimate
##     loglik.new <- logP.Ak.omega(A.new, omega)

##     p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
##     b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

##     # finalize update of A0
##     A0[b.switch,] <- A.new[b.switch,]

##     ### update Aks ###
##     deltas <- abs(Ak[,1,] - Ak[,2,]) / 2
##     lims <- A0 # this tells us how big delta can be before going outside (0,1)
##     lims <- ifelse(lims > .5, 1 - lims, lims)

##     # bring all deltas within the limits
##     deltas <- ifelse(deltas > lims, lims, deltas)

##     loglik.cur <- logP.delta.sigma(deltas, sigma, lims, dev)

##     # sample new candidate for delta
##     A1.new <- rdirichlet(dim(A0)[1], omega + totals1 + 0.5)
##     A2.new <- rdirichlet(dim(A0)[1], omega + totals2 + 0.5)

##     deltas.new <- (abs(A1.new - A0) + abs(A2.new - A0)) / 2

##     # make decision to keep new or old estimate
##     loglik.new <- logP.delta.sigma(deltas.new, sigma, lims, dev)

##     p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)

##     if(any(is.nan(p.switch))) # if we get any changes that are out of whack (too radical), reset
##     {
##         deltas.new[is.nan(p.switch),] <- 0
##         p.switch[is.nan(p.switch)] <- 1
##     }

##     b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

##     # finalize update of Ak (dependent on updated values for A and delta)
##     deltas[b.switch,] <- deltas.new[b.switch,]
##     Ak[b.switch,1,] <- A1.new[b.switch,]
##     Ak[b.switch,2,] <- A2.new[b.switch,]

##     Ak[,1,] <- A0 + ifelse(Ak[,1,] > Ak[,2,], 1, -1) * deltas
##     Ak[,2,] <- A0 + ifelse(Ak[,1,] <= Ak[,2,], 1, -1) * deltas

##     ### update sigma ###
##     sigma <- updt.sigma(sigma, deltas, lims, dev)

##     ### probability of no admixture ###
##     P.Aeq1 <- Ak # same dimensions
##     P.Aeq1[,,] <- NA

##     # this will be beta distributed for each individual population (e.g. African vs something else)
##     if(dim(P.Aeq1)[3] == 2)
##     {
##         for(l in 1:2)
##         {
##             P.Aeq1[,1,l] <- pbeta(1 - epsilon, omega[l] + totals1[,l],
##                                   omega[-l] + totals1[,-l], lower.tail = FALSE)
##             P.Aeq1[,2,l] <- pbeta(1 - epsilon, omega[l] + totals2[,l],
##                                   omega[-l] + totals2[,-l], lower.tail = FALSE)
##         }
##     }else{
##         for(l in 1:dim(P.Aeq1)[3])
##         {
##             P.Aeq1[,1,l] <- pbeta(1 - epsilon, omega[l] + totals1[,l],
##                                   sum(omega[-l]) + apply(totals1[,-l], 1, sum), lower.tail = FALSE)
##             P.Aeq1[,2,l] <- pbeta(1 - epsilon, omega[l] + totals2[,l],
##                                   sum(omega[-l]) + apply(totals2[,-l], 1, sum), lower.tail = FALSE)
##         }
##     }

##     return(list(A0 = A0, Ak = Ak, sigma = sigma, P.Aeq1 = P.Aeq1))
## }


###########################################################################################
# Update A0 and Aks, given all local ancestry estimates and distances between the markers #
###########################################################################################

updt.AX <- function(A0, AX, Ak, Ajs, omegaX, sigma, W, dev)
{
    warning('This function is depricated!')
    totals1 <- apply(Ajs[,,1,], c(1,3), sum) # summing over markers
    totals2 <- apply(Ajs[W,,2,], c(1,3), sum)
    totals <- totals1
    totals[W,] <- totals[W,] + totals2

    AX0 <- apply(AX, c(1,3), mean, na.rm = TRUE)
    A <- A0
    A[W,] <- Ak[W,1,]

    loglik.cur <- logP.AX.omega(AX0, A, omegaX)

    ### update AX0 ###
    AX0.new <- rdirichlet(dim(AX)[1], omegaX*A + totals + 0.5)

    loglik.new <- logP.AX.omega(AX0.new, A, omegaX)

    # make decision for each individual
    p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
    b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

    # make and return final selection
    AX0[b.switch,] <- AX0.new[b.switch,]

    ### update AX ###
    deltas <- abs(AX[W,1,] - AX[W,2,]) / 2
    lims <- AX0[W,] # this tells us how big delta can be before going outside (0,1)
    lims <- ifelse(lims > .5, 1 - lims, lims)

    # bring all deltas within limits
    deltas <- ifelse(deltas > lims, lims, deltas)

    loglik.cur <- logP.delta.sigma(deltas, sigma, lims, dev)

    # sample new candidate for delta
    AX1.new <- rdirichlet(length(W), omegaX*A[W,] + totals1[W,] + 0.5)
    AX2.new <- rdirichlet(length(W), omegaX*A[W,] + totals2 + 0.5)

    deltas.new <- (abs(AX1.new - AX0[W,]) + abs(AX2.new - AX0[W,])) / 2

    # make decision to keep new or old estimate
    loglik.new <- logP.delta.sigma(deltas.new, sigma, lims, dev)

    p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)

    if(any(is.nan(p.switch))) # if we get any changes that are out of whack (too radical), reset
    {
        deltas.new[is.nan(p.switch),] <- 0
        p.switch[is.nan(p.switch)] <- 1
    }

    b.switch <- as.logical(rbinom(length(p.switch), 1, p.switch))

    # finalize update of Ak (dependent on updated values for A and delta)
    deltas[b.switch,] <- deltas.new[b.switch,]
    AX[W[b.switch],1,] <- AX1.new[b.switch,]
    AX[W[b.switch],2,] <- AX2.new[b.switch,]

    AX[W,1,] <- AX0[W,] + ifelse(AX[W,1,] > AX[W,2,], 1, -1) * deltas
    AX[W,2,] <- AX0[W,] + ifelse(AX[W,1,] <= AX[W,2,], 1, -1) * deltas

    if(length(W) != dim(AX0)[1]) # just in case they are all women
        AX[-W,1,] <- AX0[-W,]

    return(AX)
}


############
# Update G #
############

## # returns an array of ancestral probabilities (two parents combined) with dimensions:
## #   i - individual
## #   j - marker
## #   o - combination of populations inherited (one from each parent)
## updt.Aj <- function(gammas) # I think this one is OK with X data, the way we have gammas.joint calculated
## {
##     n <- dim(gammas)[4]
##     n.combos <- n * (n+1) / 2

##     Aj <- array(NA, dim = c(dim(gammas)[1:2], n.combos),
##                 dimnames = list(dimnames(gammas)[[1]], dimnames(gammas)[[2]],
##                                 rep('g', n.combos)))

##     curr <- 1
##     for(l in 1:n)
##     {
##         for(l.alt in l:n)
##         {
##             if(l == l.alt){
##                 # keep the homozygous ancestral states the way they are
##                 Aj[,,curr] <- gammas[,,1,l] * gammas[,,2,l.alt]
##             }else{
##                 # add states of heterozygous states -- could go either way
##                 Aj[,,curr] <- gammas[,,1,l] * gammas[,,2,l.alt] +
##                               gammas[,,1,l.alt] * gammas[,,2,l]
##             }

##             dimnames(Aj)[[3]][curr] <- paste('g', l, l.alt, sep = '')
##             curr <- curr + 1
##         }
##     }

##     return(Aj)
## }

## updt.G0 <- function(G)
## {
##     n <- dim(G)[4]
##     n.combos <- n * (n+1) / 2

##     G0 <- array(NA, dim = c(dim(G)[1:2], n.combos),
##                 dimnames = list(dimnames(G)[[1]], dimnames(G)[[2]], rep('g', n.combos)))

##     curr <- 1
##     for(l in 1:n)
##     {
##         for(l.alt in l:n)
##         {
##             G0[,,curr] <- (G[,,1,l] == 1 & G[,,2,l.alt] == 1) |
##                           (G[,,2,l] == 1 & G[,,1,l.alt] == 1)

##             dimnames(G0)[[3]][curr] <- paste('g', l, l.alt, sep = '')
##             curr <- curr + 1
##         }
##     }

##     return(G0)
## }

updt.G <- function(G, G.new, b.switch)
{
    # if G is NULL, we're on the first round -> keep it no matter how bad
    if(is.null(G))
    {
        return(G.new)
    }

    for(i in 1:dim(G)[1])
    {
        for(k in 1:dim(G)[3])
        {
            index <- i + (k - 1) * dim(G)[3]
            if(b.switch[index] == 1)
                G[i,,k,] <- G.new[i,,k,]
        }
    }

    return(G)
}
