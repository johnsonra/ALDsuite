# admixture.par.R
# code for the HMM, which can be easily loaded for parallel execution
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created April 17, 2013
# Last Modified March 19, 2014


admixture.par <- function()
{
    return(quote({

    ####################
    # Calculate Gammas #
    ####################

    point <- 0
    if(!is.null(haps))
    {
        gammas <- P.gammas(haps, P, Pm.prior, lambda, lambdaX, d, chr,
                           A0, Ak, AX, gender, sex.chr, haps = TRUE, dev)
    }else{
        gammas <- P.gammas(geno, P, Pm.prior, lambda, lambdaX, d, chr, A0, Ak, AX, gender, sex.chr,
                           haps = FALSE, dev)
    }


    ############
    # Sample G #
    ############

    if(!exists("G"))
    {
        loglik.cur <- gammas[,1,,1]
        loglik.cur[] <- -Inf
        G <- NULL
    }

    point <- 1

    # Metropolis Hastings sample of G
    for(cyc in 1:every)
    {
        # sample ancestral states based on likelihoods of gammas
        G.new <- sample.G(gammas, lambda, d)

        # Calculate likelihood of G conditional on the HMM
        if(is.null(haps))
        {
            loglik.new <- logP.G.lambda.d(G.new, Ak, lambda, d, phased = FALSE)

            p.switch <- sapply(exp(loglik.new - loglik.cur), min, 1)
            b.switch <- sapply(p.switch, rbinom, n = 1, size = 1)
        }else{
            loglik.new <- logP.G.lambda.d(G.new, Ak, lambda, d, phased = TRUE)

            p.switch <- apply(exp(loglik.new - loglik.cur), 1:2, min, 1)
            b.switch <- apply(p.switch, 1:2, rbinom, n = 1, size = 1)
        }

        # Keep the updates or reject them
        loglik.cur[as.logical(b.switch)] <- loglik.new[as.logical(b.switch)]
        G <- updt.G(G, G.new, b.switch)

        # if even or the last cycle, sample
        if(cyc %% 2 == 0 | cyc == every)
        {
            if(!exists("G.tmp"))
            {
                G.tmp <- G
                n <- 1
            }else{
                G.tmp <- G.tmp + G
                n <- n + 1
            }
        }
    }

    # update gammas
    if(nburn + niter < 1)
    {
        gammas <- G.tmp / n
        rm(G.tmp)
        
        # take maximum likelihood estimate of G for update of parameters? Less noisy...
        G <- ML.G(gammas)
    }


    point <- 4
    ###############
    # Update Step #
    ###############

    #################### Primary Parameters ####################

    ######### A #########
    tmp <- updt.A(A0, Ak, G[, chr != sex.chr,,], omega, sigmaA, epsilon = 0.01, !is.null(haps), dev)
    A0 <- tmp$A0
    Ak <- tmp$Ak
    sigmaA <- tmp$sigma
    P.Aeq1 <- tmp$P.Aeq1

    ## point <- 5
    ## ######### AX #########
    ## if(any(chr == sex.chr))
    ## {
    ##     AX <- updt.AX(A0, AX, Ak, G[, chr == sex.chr,,], omegaX, sigmaA, W, dev)

    ##     # only swap these out during the burn-in phase
    ##     if(burnin)
    ##     {
    ##         # check how well X and autosome ancestries match up in males
    ##         mom <- apply(abs(Ak[M, 1,] - AX[M, 1,]), 1, sum)
    ##         dad <- apply(abs(Ak[M, 2,] - AX[M, 1,]), 1, sum)

    ##         swap <- mom > dad # more of a deviation on mom than in dad...

    ##         tmpA <- Ak[M[swap],1,]
    ##         tmplambda <- lambda[M[swap],1]

    ##         Ak[M[swap],1,] <- Ak[M[swap],2,]
    ##         lambda[M[swap],1] <- lambda[M[swap],2]
    ##     }

    ##     AX0 <- apply(AX, c(1,3), mean, na.rm = TRUE)
    ## }else{
    ##     AX[,,] <- NA
    ## }

    point <- 6
    ######### lambda #########

    lambda <- updt.lambda(G[,chr != sex.chr,,, drop = FALSE], Ak, lambda, Pm.prior, d[chr != sex.chr],
                          alpha[1], alpha[2], chr[chr != sex.chr], phased = !is.null(haps), dev)

    ## point <- 7
    ## ######### lambdaX #########
    ## if(any(chr == sex.chr))
    ## {
    ##     tmp <- updt.lambda(G[,chr == sex.chr,,], G0[,chr == sex.chr,], AX0, AX,
    ##                        lambdaX, sigmaL, d[chr == sex.chr], alpha[1], alpha[2],
    ##                        chr[chr == sex.chr], P.Aeq1, dev, W)
    ##     lambdaX <- tmp$lambda
    ## }else{
    ##     lambdaX[,] <- NA
    ## }

    point <- 8
    ######### P #########
    # always update this during the first few cycles
    # after that, update with frequency dependent upon tau
    if(supercyc <= 5 | rbinom(1, 1, 1 / max(1, min(sqrt(tau)))))
    {
        if(!is.null(haps))
        {
            tmp <- updt.P.betas(P, Pm, Pm.prior, lambda, tau, G, haps, phased = TRUE, dev)
        }else{
            tmp <- updt.P.betas(P, Pm, Pm.prior, lambda, tau, G, geno, phased = FALSE, dev)
        }

        P <- tmp$P
        Pm.prior <- tmp$Pm.prior

        point <- 9
    #################### Hyper Parameters ####################

        # don't update every time -- more often when tau is small
        ######### tau and Pm #########
        if(supercyc == 1) # initialize this the first time through... no need to repeat
            to.update <- updt.tau.Pm.prep(tau, Pm.prior)

        tmp <- updt.tau.Pm(tau, Pm, Pm.counts, Pm.prior, P, to.update, dev = dev)
        tau <- tmp$tau
        Pm <- tmp$Pm
        Pm.prior <- tmp$Pm.prior
    }

    point <- 10
    ######### alpha and beta #########
    tmp <- updt.alpha.beta(alpha[1], alpha[2], apply(lambda, 1, mean))
    alpha <- c(tmp$alpha, tmp$beta)

    ## point <- 11
    ## if(any(chr %in% sex.chr))
    ## {
    ##     tmp <- updt.alpha.beta(alphaX[1], alphaX[2], apply(lambda, 1, mean, na.rm = TRUE))
    ##     alphaX <- c(tmp$alpha, tmp$beta)
    ## }else{
    ##     alphaX <- rep(NA, 2)
    ## }

    point <- 12
    ######### omega #########
    omega <- updt.omega(omega, A0)

    ## point <- 13
    ## if(any(chr == sex.chr))
    ## {
    ##     omegaX <- updt.omegaX(omegaX, Ak, AX, W)
    ## }else{
    ##     omegaX <- NA
    ## }

    point <- 14
    if(!burnin)
    {
        final$A0 <- final$A0 + A0
        final$Ak <- final$Ak + Ak
        final$gammas <- final$gammas + G
        final$lambda <- final$lambda + lambda
        final$lambdaX <- final$lambdaX + lambdaX
        final$P <- final$P + P
        final$Pm <- final$Pm + Pm
        final$omega <- final$omega + omega
        final$omegaX <- final$omegaX + omegaX
        final$alpha <- final$alpha + alpha
        final$alphaX <- final$alphaX + alphaX
        final$tau <- final$tau + tau
        final$sigmaA <- final$sigmaA + sigmaA
        final$sigmaL <- final$sigmaL + sigmaL
    }
  }))
}
