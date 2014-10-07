# admixture.R
# Infer ancestry estimates for a data set
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created September 10, 2012
# Last Modified October 7, 2014

# see documentation for details on these arguments
admixture <- function(Pm.prior, haps = NULL, geno = NULL, gender = NULL, chr, pos, burn = 100,
                      iter = 200, every = 5, indiv.id = NULL, marker.id = NULL, pop.id = NULL,
                      lambda = 9, tau = 300, omega = NULL, rand.seed = NULL,
                      cores = detectCores(), cl = NULL, dev = FALSE, verbose = TRUE,
                      debug = FALSE, sex.chr = 23, indiv.geno.check = 0.98,
                      marker.geno.check = .98, hwe.thresh = 1e-4, bad.indiv = NULL,
                      male.het.X = 1, fast = FALSE, A0 = NULL, Ak = NULL, run.checks = TRUE)
{

######### haplo/geno-type formatting #########

    if(!is.null(haps) & is.null(geno))
    {
        tmp <- array(0, dim = c(dim(haps)[1] / 2, dim(haps)[2], 2),
                     dimnames = list(dimnames(haps)[[1]][1:(dim(haps)[1]/2) * 2 - 1],
                                     dimnames(haps)[[2]], c('Mother', 'Father')))

        tmp[,,1] <- haps[1:dim(tmp)[1] * 2 - 1,]
        tmp[,,2] <- haps[1:dim(tmp)[1] * 2,]
        haps <- tmp

        geno <- apply(haps, 1:2, sum)
    }

    # check that dimnames, and IDs are consistent
    if(!is.null(haps))
    {
        if(is.null(dimnames(haps)[[1]]))
        {
            if(!is.null(indiv.id))
            {
                dimnames(haps)[[1]] <- indiv.id
            }else{
                warning("Neither indiv.id nor dimnames of haps defined.")
                dimnames(haps)[[1]] <- 1:dim(haps)[1]
            }
        }
    }


######### QC checks on data formatting and genetic data #########
    if(run.checks)
        ald.qc(Pm.prior, haps, geno, gender, chr, pos, burn,
               iter, every, indiv.id, marker.id, pop.id,
               lambda, tau, omega, rand.seed,
               cores, cl, dev, verbose,
               debug, sex.chr, indiv.geno.check,
               marker.geno.check, hwe.thresh, bad.indiv,
               bad.marker, male.het.X, fast, A0, Ak)


#########
# Setup #
#########


######### Gender #########
    W <- which(gender == 'F')
    M <- which(gender == 'M')


######### Distances #########
    # check this to be sure that we have the correct distance measure!!!
    if(max(pos) > 1000)
    {
        warning('Physical position given rather than genetic position? Assuming 1 cM ~ 1 Mb.')
        divisor <- 1e6
    }else{
        if(max(pos) > 10)
        {
            divisor <- 1 # assume possition is in cM
        }else{
            warning('Assuming genetic is position in Morgans.')
            divisor <- 0.1 # probably in morgans, otherwise
        }
    }

    # distance to next marker (in Centimorgans)
    d <- c(NA, pos[-1] - pos[-length(pos)]) / divisor

    # start each chromosome wth an NA
    d <- ifelse(c(TRUE, chr[-1] != chr[-length(chr)]), NA, d)

    # make sure pos is ordered properly
    ##### Add this to the checks...not here #####
    # -- there may be some slight rounding error for very close markers, but there shouldn't be
    if(any(pos[-1] < pos[-length(pos)]))
    {
        diffs <- pos[-1] - pos[-length(pos)]
        bad <- which(diffs < 0)

        ## if(any(diffs[bad] < -0.00000001)) # this should be fixed now...
        ##     stop("Some positions appear not to be ordered correctly.")

        d <- abs(d)
    }


######### starting estimate of ancestral allele frequencies #########
    P <- matrix(NA, nrow = length(Pm.prior), ncol = length(Pm.prior[[1]]$freq),
                dimnames = list(marker.id, pop.id))

    Pm.counts <- array(dim = c(dim(P), 2),
                       dimnames = list(dimnames(P)[[1]], dimnames(P)[[2]],
                                       c('ref', 'vart')))

    for(j in 1:dim(P)[1])
    {
        for(k in 1:dim(P)[2])
        {
            Pm.counts[j,k,'vart'] <- Pm.prior[[j]]$freq[k] * Pm.prior[[j]]$n[k] + 0.5
            Pm.counts[j,k,'ref'] <- (1 - Pm.prior[[j]]$freq[k]) * Pm.prior[[j]]$n[k] + 0.5

            if(Pm.counts[j,k,'vart'] + Pm.counts[j,k,'ref'] == 1) # no prior data -> switch to Unif(0,1)
            {
                Pm.counts[j,k,'vart'] <- 1
                Pm.counts[j,k,'ref'] <- 1
            }
        }
    }

    P <- Pm.counts[,,'vart'] / apply(Pm.counts, 1:2, sum)
    Pm <- P

    # keep this for future use -- P will change
    P.ori <- P

######### starting estimate of ancestry ######### !!!!!!!! this is far from eligant, but it is fast
#                                                          - make it better!
    ## if(is.null(A0))
    ## {
        A0 <- matrix(0, nrow = dim(geno)[1], ncol = dim(P)[2])

        for(j in 1:dim(P)[1])
        {
        # these are the ancestral combinations
        combos <- matrix(nrow = 0, ncol = 2)
        for(l in 1:dim(P)[2])
            combos <- rbind(combos, cbind(l, l:dim(P)[2]))

        # this is what we multiply by (1 if they are the same, two if they are different
        # -- eg could be {1, 2} or {2, 1})
            combos <- cbind(combos, ifelse(combos[,1] == combos[,2], 1, 2))

            # this is where we'll collect all the p-values for each posibility
            ps <- matrix(0, nrow = dim(geno)[1], ncol = dim(combos)[1])
            for(g in 1:dim(combos)[1])
                ps[,g] <- P.a.gamma.gamma(geno[,j], combos[g,1], combos[g,2], P[j,]) * combos[g,3]

            # pick the column with the maximum probability
            maxs <- apply(ps, 1, which.max)

            # take care of NAs -- replace with 0 -- thus they won't be counted below
            if(class(maxs) == 'list')
            {
                maxs <- sapply(maxs, function(x)
                           {
                               if(length(x) == 0)
                                   return(0)
                               return(x)
                           })
            }

            # count up how many of each ancestral population we have
            for(l in 1:dim(P)[2])
            {
                A0[combos[maxs,1] == l,l] <- A0[combos[maxs,1] == l,l] + 1
                A0[combos[maxs,2] == l,l] <- A0[combos[maxs,2] == l,l] + 1
            }
        }

        # count up how many total we have in each row
        tmp <- apply(A0, 1, sum)
        A0 <- A0 / tmp # individual's global ancestry estimates
        dimnames(A0) <- list(indiv.id, pop.id)
    ## }

    ## if(is.null(Ak))
    ## {
        Ak <- array(dim = c(dim(A0)[1], 2, dim(P)[2]),
                    dimnames = list(indiv.id, c('Mother', 'Father'), pop.id))
        Ak[,1,] <- A0 # global ancestry estimates for parents
        Ak[,2,] <- A0
    ## }

    # X chromosomes
    AX <- Ak
    AX[M, 'Father',] <- NA

######### starting estimate of tau #########
    tau <- rep(tau, dim(P)[2])[1:dim(P)[2]] # in case they gave us a tau value for each populaton

########## starting estimate of lambda #########
    lambda <- matrix(lambda, nrow = dim(geno)[1], ncol = 2,
                     dimnames = list(indiv.id, c('Mother', 'Father'))) # one for each parent

    # make lambdaX NA for males
    lambdaX <- lambda
    lambdaX[M, 'Father'] <- NA

######### starting estimates for beta and sigma #########

    maxbetas <- 1

    for(j in 1:length(Pm.prior))
    {
        for(k in 1:length(Pm.prior[[j]]$model))
        {
            # we use the upstream markers to infer genotypes in forward chain
            if(length(Pm.prior[[j]]$model[[k]]$betas) > 0)
            {
                Pm.prior[[j]]$model[[k]]$betas.p <- Pm.prior[[j]]$model[[k]]$betas
                Pm.prior[[j]]$model[[k]]$betas.prior <- Pm.prior[[j]]$model[[k]]$betas

                Pm.prior[[j]]$model[[k]]$sigma <- diag(Pm.prior[[j]]$model[[k]]$vcv)

                if(any(is.na(Pm.prior[[j]]$model[[k]]$sigma)))
                    stop("NaNs detected in Pm.prior[[...]]model[[...]]$vcv")

                if(length(Pm.prior[[j]]$model[[k]]$betas) > maxbetas)
                    maxbetas <- length(Pm.prior[[j]]$model[[k]]$betas)
            }
        }
    }


########## hyper parameters ######### !!!!!!! these could be tuned some!
    alpha <- c(lambda[1,1], 1)

    alphaX <- alpha

    if(is.null(omega))
        omega <- apply(A0, 2, mean) * 10 # this will be forced to be >= 1 later

    omegaX <- 1

    sigmaA <- rep(.2, length(omega))
    names(sigmaA) <- names(omega)

    sigmaL <- 3

#################### values to save ####################

    # generate names for Aj
    n <- dim(P)[2]

    combos <- character(n * (n+1) / 2)

    curr <- 1
    for(i in 1:n)
    {
        for(j in i:n)
        {
            combos[curr] <- paste('g', i, j, sep = '')
            curr <- curr + 1
        }
    }

    # this is the number of burn-in samples / follow-on iterations we will end up with
    nburn <- ceiling(burn / cores)

    niter <- ceiling(iter / cores)

    ## if(verbose)
    ## {
    ##     burn.lik <- list(Ak = numeric(), AX = numeric(), Aj = numeric(), lambda = numeric(),
    ##                      lambdaX = numeric(), P = numeric(), Pm = numeric())
    ##     iter.lik <- list(Ak = numeric(), AX = numeric(), Aj = numeric(), lambda = numeric(),
    ##                      lambdaX = numeric(), P = numeric(), Pm = numeric())
    ## }else{
        burn.lik <- NULL
        iter.lik <- NULL
    ## }

    if(debug & cores == 1)
    {
        tmp <- require(ncdf4)

        if(!tmp)
        {
            warning("ncdf package required to run in debug mode.")
            debug <- FALSE
        }else{
            # set up dimensions
            dimIndivs <- ncdim_def("Indivs", "count", 1:dim(haps)[1])
            dimMarkers <- ncdim_def("Markers", "count", 1:dim(haps)[2])
            dimParents <- ncdim_def("Parents", "count", 1:2)
            dimPops <- ncdim_def("Pops", "count", 1:dim(Pm.counts)[3])
            dimIter <- ncdim_def("Iter", "count", 1:(burn + iter))
            dimHyper2 <- ncdim_def("Hyper2", "count", 1:2) # for hyper parameters with two values (like omega)
            dimBetas <- ncdim_def("Betas", "count", 1:maxbetas)

            # declare variables
            debug.vars <- list(ncvar_def("A", "ancestry", list(dimIndivs, dimPops, dimIter), -1),
                               ncvar_def("Ak", "ancestry", list(dimIndivs, dimParents, dimPops, dimIter), -1),
                               ncvar_def("AX", "ancestry", list(dimIndivs, dimParents, dimPops, dimIter), -1),
                               ncvar_def("lambda", "generation", list(dimIndivs, dimParents, dimIter), -1),
                               ncvar_def("lambdaX", "generation", list(dimIndivs, dimParents, dimIter), -1),
                               ncvar_def("P", "frequency", list(dimMarkers, dimPops, dimIter), -1),
                               ncvar_def("Pm", "frequency", list(dimMarkers, dimPops, dimIter), -1),
                               ncvar_def("omega", "hyper", list(dimHyper2, dimIter), -1),
                               ncvar_def("omegaX", "hyper", list(dimIter), -1),
                               ncvar_def("alpha", "hyper", list(dimHyper2, dimIter), -1), # not really parents, but has two anyway
                               ncvar_def("alphaX", "hyper", list(dimParents, dimIter), -1),
                               ncvar_def("tau", "hyper", list(dimHyper2, dimIter), -1),
                               ncvar_def("sigmaA", "hyper", list(dimPops, dimIter), -1),
                               ncvar_def("sigmaL", "hyper", list(dimIter), -1),
                               ncvar_def("BetasUp", "hyper", list(dimBetas, dimMarkers, dimPops, dimIter), -1),
                               ncvar_def("BetasDown", "hyper", list(dimBetas, dimMarkers, dimPops, dimIter), -1))

            # create and open file
            nc.debug <- nc_create("mald.debug.nc", debug.vars)

            nc.debug <- nc_close(nc.debug)
        }
    }

    if(debug & cores > 1)
    {
        debug <- FALSE
        warning("debug mode only supported when cores == 1\nturning debug mode off")
    }


##################
# Set up Cluster #
##################

    hmm <- admixture.par()

    burnin <- TRUE # this is what we start with

    if(cores > 1)
    {
        if(is.null(cl))
        {
            cl <- makeCluster(cores)
            stop.cl <- TRUE
        }else{
            stop.cl <- FALSE
        }

        clusterSetRNGStream(cl, rand.seed)

        invisible(clusterEvalQ(cl, library(ALDsuite)))

        clusterExport(cl, list('hmm', 'haps', 'geno', 'P', 'lambda', 'lambdaX', 'd', 'chr', 'A0', 'Ak', 'AX',
                               'gender', 'sex.chr', 'dev', 'omega', 'omegaX', 'W', 'M', 'alpha', 'alphaX',
                               'Pm', 'Pm.counts', 'Pm.prior', 'every', 'tau', 'verbose', 'cores',
                               'burnin', 'burn.lik', 'iter.lik', 'sigmaA', 'sigmaL', 'debug'),
                      envir = environment())

        # get system pids for each subprocess
        ids <- unlist(clusterEvalQ(cl, id <- Sys.getpid()))
        clusterExport(cl, list('ids'))

        # now relabel them from 1:cores
        clusterEvalQ(cl, id <- which(id == ids))

    }else{
        stop.cl <- FALSE
    }


##############################################################################
#################################### MCMC ####################################
##############################################################################

    if(nburn > 0)
    {
        iter.every <- every

        if(cores > 1)
            invisible(clusterEvalQ(cl, every <- 1))

        for(supercyc in 1:(nburn / every))
        {
            if(cores > 1)
            {
                clusterExport(cl, list('supercyc'), envir = environment())

                if(verbose) # report on progress
                    invisible(sapply(clusterEvalQ(cl, (supercyc - 1) * cores + id), cat, '\n'))

                invisible(clusterEvalQ(cl, eval(hmm)))
            }else{

                if(verbose) # report on progress
                    cat(supercyc, '\n')

                eval(hmm, envir = environment())

                if(debug)
                {
                    nc.debug <- nc_open('mald.debug.nc', write = TRUE)

                    # write values
                    ncvar_put(nc.debug, "A", A0, start = c(1, 1, supercyc), count = c(dim(A0), 1))
                    ncvar_put(nc.debug, "Ak", Ak, start = c(1, 1, 1, supercyc), count = c(dim(Ak), 1))
                    ncvar_put(nc.debug, "AX", AX, start = c(1, 1, 1, supercyc), count = c(dim(AX), 1))
                    ncvar_put(nc.debug, "lambda", lambda, start = c(1, 1, supercyc), count = c(dim(lambda), 1))
                    ncvar_put(nc.debug, "lambdaX", lambdaX, start = c(1, 1, supercyc),
                              count = c(dim(lambdaX), 1))
                    ncvar_put(nc.debug, "P", P, start = c(1, 1, supercyc), count = c(dim(P), 1))
                    ncvar_put(nc.debug, "Pm", Pm, start = c(1, 1, supercyc), count = c(dim(Pm), 1))
                    ncvar_put(nc.debug, "omega", omega, start = c(1, supercyc), count = c(length(omega), 1))
                    ncvar_put(nc.debug, "omegaX", omegaX, start = supercyc, count = 1)
                    ncvar_put(nc.debug, "alpha", alpha, start = c(1, supercyc), count = c(2, 1))
                    ncvar_put(nc.debug, "alphaX", alphaX, start = c(1, supercyc), count = c(2, 1))
                    ncvar_put(nc.debug, "tau", tau, start = c(1, supercyc), count = c(length(tau), 1))
                    ncvar_put(nc.debug, "sigmaA", sigmaA, start = c(1, supercyc), count = c(length(sigmaA), 1))
                    ncvar_put(nc.debug, "sigmaL", sigmaL, start = supercyc, count = 1)

                    for(j in 1:length(Pm.prior))
                    {
                        for(k in 1:length(Pm.prior[[1]]$upstream))
                        {
                            tmp <- rep(0, maxbetas)
                            b <- Pm.prior[[j]]$upstream[[k]]$betas
                            if(length(b) > 0)
                                tmp[1:length(b)] <- b
                            ncvar_put(nc.debug, "BetasUp", tmp, start = c(1, j, k, supercyc),
                                      count = c(maxbetas, 1, 1, 1))

                            tmp <- rep(0, maxbetas)
                            b <- Pm.prior[[j]]$downstream[[k]]$betas
                            if(length(b) > 0)
                                tmp[1:length(b)] <- b
                            ncvar_put(nc.debug, "BetasDown", tmp, start = c(1, j, k, supercyc),
                                      count = c(maxbetas, 1, 1, 1))
                        }
                    }
                    nc.debug <- nc_close(nc.debug)
                }
            }

            # this is for the *next burn-in* sample (i.e. skip the last cycle)
            if(cores > 1 & supercyc != max(1:(nburn / every)))
            {
                ##### get local state for each chain #####

                # global ancestry
                Ak.r <- clusterEvalQ(cl, Ak)
                AX.r <- clusterEvalQ(cl, AX)

                omega.r <- clusterEvalQ(cl, omega)
                omegaX.r <- mean(unlist(clusterEvalQ(cl, omegaX)))

                # mean generations since admixture event
                lambda.r <- clusterEvalQ(cl, lambda)
                lambdaX.r <- clusterEvalQ(cl, lambdaX)

                alpha.r <- clusterEvalQ(cl, alpha)
                alphaX.r <- clusterEvalQ(cl, alphaX)

                # ancestry specific allele frequencies
                P.r <- clusterEvalQ(cl, P)
                Pm.r <- clusterEvalQ(cl, Pm)
                Pm.prior.all <- clusterEvalQ(cl, Pm.prior) ############### this could be optimized a bit...don't need all of that info!

                tau.r <- clusterEvalQ(cl, tau)

                ##### calculate remote proposal #####
                Pm.prior.r <- Pm.prior.all[[1]]

                for(s in 2:cores)
                {
                    Ak.r[[1]] <- Ak.r[[1]] + Ak.r[[s]]
                    AX.r[[1]] <- AX.r[[1]] + Ak.r[[s]]
                    omega.r[[1]] <- omega.r[[1]] + omega.r[[s]]
                    alpha.r[[1]] <- alpha.r[[1]] + alpha.r[[s]]
                    alphaX.r[[1]] <- alphaX.r[[1]] + alphaX.r[[s]]
                    lambda.r[[1]] <- lambda.r[[1]] + lambda.r[[s]]
                    lambdaX.r[[1]] <- lambdaX.r[[1]] + lambdaX.r[[s]]
                    P.r[[1]] <- P.r[[1]] + P.r[[s]]
                    Pm.r[[1]] <- Pm.r[[1]] + Pm.r[[s]]
                    tau.r[[1]] <- tau.r[[1]] + tau.r[[s]]

                    for(j in 1:length(Pm.prior))
                    {
                        for(k in 1:length(Pm.prior[[j]]$model))
                        {
                            if(length(Pm.prior[[j]]$model[[k]]$linked) > 0)
                            {
                                Pm.prior.r[[j]]$model[[k]]$betas <- Pm.prior.r[[j]]$model[[k]]$betas +
                                                                    Pm.prior.all[[s]][[j]]$model[[k]]$betas
                                Pm.prior.r[[j]]$model[[k]]$vcv <- Pm.prior.r[[j]]$model[[k]]$vcv +
                                                                  Pm.prior.all[[s]][[j]]$model[[k]]$vcv
                                Pm.prior.r[[j]]$model[[k]]$hessian <- Pm.prior.r[[j]]$model[[k]]$hessian +
                                                                      Pm.prior.all[[s]][[j]]$model[[k]]$hessian
                            }
                        }
                    }
                }

                Ak.r <- Ak.r[[1]] / cores
                AX.r <- AX.r[[1]] / cores
                omega.r <- omega.r[[1]] / cores
                lambda.r <- lambda.r[[1]] / cores
                lambdaX.r <- lambdaX.r[[1]] / cores
                alpha.r <- alpha.r[[1]] / cores
                alphaX.r <- alphaX.r[[1]] / cores
                P.r <- P.r[[1]] / cores
                Pm.r <- Pm.r[[1]] / cores
                tau.r <- tau.r[[1]] / cores

                for(j in 1:length(Pm.prior))
                {
                    for(k in 1:length(Pm.prior[[j]]$model))
                    {
                        if(length(Pm.prior[[j]]$model[[k]]$linked) > 0)
                        {
                            Pm.prior.r[[j]]$model[[k]]$betas <- Pm.prior.r[[j]]$model[[k]]$betas / cores
                            Pm.prior.r[[j]]$model[[k]]$vcv <- Pm.prior.r[[j]]$model[[k]]$vcv / cores
                            Pm.prior.r[[j]]$model[[k]]$hessian <- Pm.prior.r[[j]]$model[[k]]$hessian / cores
                        }
                    }
                }

                ##### update local proposal states #####
                mixing <- supercyc / max(1:(nburn/every)) # This is for the next cycle -- first time is close to 0

                clusterExport(cl, list('Ak.r', 'AX.r', 'omega.r', 'omegaX.r', 'lambda.r',
                                       'lambdaX.r', 'alpha.r', 'alphaX.r',
                                       'P.r', 'Pm.r', 'Pm.prior.r', 'tau.r', 'mixing'),
                              envir = environment())

                invisible(clusterEvalQ(cl,
                             {
                                 Ak <- mixing*Ak + (1 - mixing)*Ak.r
                                 AX <- mixing*AX + (1 - mixing)*AX.r

                                 omega <- mixing*omega + (1 - mixing)*omega.r
                                 omegaX <- mixing*omegaX + (1 - mixing)*omegaX.r

                                 lambda <- mixing*lambda + (1 - mixing)*lambda.r
                                 lambda.X <- mixing*lambdaX + (1 - mixing)*lambdaX.r

                                 alpha <- mixing*alpha + (1 - mixing)*alpha.r
                                 alphaX <- mixing*alphaX + (1 - mixing)*alphaX.r

                                 P <- mixing*P + (1 - mixing)*P.r
                                 Pm <- mixing*Pm + (1 - mixing)*Pm.r

                                 tau <- mixing*tau + (1 - mixing)*tau.r

                                 for(j in 1:length(Pm.prior))
                                 {
                                     for(k in 1:length(Pm.prior[[j]]$model))
                                     {
                                         if(length(Pm.prior[[j]]$model[[k]]$linked) > 0)
                                         {
                                             Pm.prior[[j]]$model[[k]]$betas <- mixing*Pm.prior[[j]]$model[[k]]$betas +
                                                                       (1 - mixing)*Pm.prior.r[[j]]$model[[k]]$betas
                                             Pm.prior[[j]]$model[[k]]$vcv <- mixing*Pm.prior[[j]]$model[[k]]$vcv +
                                                                     (1 - mixing)*Pm.prior.r[[j]]$model[[k]]$vcv
                                             Pm.prior[[j]]$model[[k]]$hessian <- mixing*Pm.prior[[j]]$model[[k]]$hessian +
                                                                         (1 - mixing)*Pm.prior.r[[j]]$model[[k]]$hessian
                                         }
                                     }
                                 }
                             }))
            }
        }
    }

    if(niter > 0)
    {
        if(exists('iter.every')) # from burnin
            every <- iter.every
        burnin <- FALSE
        final <- list(A0 = 0, Ak = 0, gammas = 0, lambda = 0, lambdaX = 0, P = 0, Pm = 0,
                      P.ori = P.ori, omega = 0, omegaX = 0, alpha = 0, beta = 0,
                      alphaX = 0, tau = 0, niter = 0, seed = rand.seed,
                      sigmaA = 0, sigmaL = 0, cores = cores)

        if(cores > 1)
        {
            # set up each chain
            invisible(clusterExport(cl, list('final', 'niter', 'burnin', 'every'), envir = environment()))

            # do all samples for each chain
            if(fast) # may not be reproducible, could be a lot faster
            {
                invisible(clusterApplyLB(cl, 1:(iter / every), function(supercyc) eval(hmm, envir = globalenv())))
            }else{

                # progress reports not yet working for parallel version
                invisible(clusterEvalQ(cl, for(supercyc in 1:(niter / every)) eval(hmm)))
            }

            # get results from each chain
            chains <- clusterEvalQ(cl, final)

            # get weighted results
            for(i in 1:cores)
            {
                final$A0 <- final$A0 + chains[[i]]$A0
                final$Ak <- final$Ak + chains[[i]]$Ak
                final$gammas <- final$gammas + chains[[i]]$gammas
                final$lambda <- final$lambda + chains[[i]]$lambda
                final$lambdaX <- final$lambdaX + chains[[i]]$lambdaX
                final$P <- final$P + chains[[i]]$P
                final$Pm <- final$Pm + chains[[i]]$Pm
                final$omega <- final$omega + chains[[i]]$omega
                final$omegaX <- final$omegaX + chains[[i]]$omegaX
                final$alpha <- final$alpha + chains[[i]]$alpha
                final$alphaX <- final$alphaX + chains[[i]]$alphaX
                final$tau <- final$tau + chains[[i]]$tau
                final$sigmaA <- final$sigmaA + chains[[i]]$sigmaA
                final$sigmaL <- final$sigmaL + chains[[i]]$sigmaL
            }

        }else{
            # running things serially here
            for(supercyc in 1:(niter / every))
            {
                if(verbose) # report on progress
                    cat(supercyc, '\n')

                eval(hmm, envir = environment())

                if(debug)
                {
                    nc.debug <- nc_open('mald.debug.nc', write = TRUE)

                    # write values
                    ncvar_put(nc.debug, "A", A0, start = c(1, 1, supercyc), count = c(dim(A0), 1))
                    ncvar_put(nc.debug, "Ak", Ak, start = c(1, 1, 1, supercyc), count = c(dim(Ak), 1))
                    ncvar_put(nc.debug, "AX", AX, start = c(1, 1, 1, supercyc), count = c(dim(AX), 1))
                    ncvar_put(nc.debug, "lambda", lambda, start = c(1, 1, supercyc), count = c(dim(lambda), 1))
                    ncvar_put(nc.debug, "lambdaX", lambdaX, start = c(1, 1, supercyc),
                              count = c(dim(lambdaX), 1))
                    ncvar_put(nc.debug, "P", P, start = c(1, 1, supercyc), count = c(dim(P), 1))
                    ncvar_put(nc.debug, "Pm", Pm, start = c(1, 1, supercyc), count = c(dim(Pm), 1))
                    ncvar_put(nc.debug, "omega", omega, start = c(1, supercyc), count = c(length(omega), 1))
                    ncvar_put(nc.debug, "omegaX", omegaX, start = supercyc, count = 1)
                    ncvar_put(nc.debug, "alpha", alpha, start = c(1, supercyc), count = c(2, 1))
                    ncvar_put(nc.debug, "alphaX", alphaX, start = c(1, supercyc), count = c(2, 1))
                    ncvar_put(nc.debug, "tau", tau, start = c(1, supercyc), count = c(length(tau), 1))
                    ncvar_put(nc.debug, "sigmaA", sigmaA, start = c(1, supercyc), count = c(length(sigmaA), 1))
                    ncvar_put(nc.debug, "sigmaL", sigmaL, start = supercyc, count = 1)

                    for(j in 1:length(Pm.prior))
                    {
                        for(k in 1:length(Pm.prior[[1]]$upstream))
                        {
                            tmp <- rep(0, maxbetas)
                            b <- Pm.prior[[j]]$upstream[[k]]$betas
                            if(length(b) > 0)
                                tmp[1:length(b)] <- b
                            ncvar_put(nc.debug, "BetasUp", tmp, start = c(1, j, k, supercyc),
                                      count = c(maxbetas, 1, 1, 1))

                            tmp <- rep(0, maxbetas)
                            b <- Pm.prior[[j]]$downstream[[k]]$betas
                            if(length(b) > 0)
                                tmp[1:length(b)] <- b
                            ncvar_put(nc.debug, "BetasDown", tmp, start = c(1, j, k, supercyc),
                                      count = c(maxbetas, 1, 1, 1))
                        }
                    }
                    nc.debug <- nc_close(nc.debug)
                }
            }
        }

        # divisor (this could be different if fast == TRUE)
        if(fast)
        {
            divisor <- floor(iter / every)
        }else{
            divisor <- cores * floor(niter / every)
        }

        final$A0 <- final$A0 / divisor
        final$Ak <- final$Ak / divisor
        final$gammas <- final$gammas / divisor
        final$lambda <- final$lambda / divisor
        final$lambdaX <- final$lambdaX / divisor
        final$P <- final$P / divisor
        final$Pm <- final$Pm / divisor
        final$omega <- final$omega / divisor
        final$omegaX <- final$omegaX / divisor
        final$alpha <- final$alpha / divisor
        final$beta <- final$beta / divisor
        final$alphaX <- final$alphaX / divisor
        final$tau <- final$tau / divisor
    }else{
        if(nburn == 0)
            gammas <- NULL

        final <- list(list(A0 = A0,
                           Ak = Ak,
                           gammas = gammas,
                           lambda = lambda,
                           lambdaX = lambdaX,
                           P = P,
                           Pm = Pm,
                           P.ori = P.ori,
                           omega = omega,
                           omegaX = omegaX,
                           alpha = alpha,
                           beta = beta,
                           alphaX = alphaX,
                           tau = tau,
                           seed = rand.seed,
                           cores = cores))
    }

##############################################################################
################################## clean up ##################################
##############################################################################

    # stop cluster if this function started it up
    if(stop.cl)
        stopCluster(cl)

    # return results
    if(!verbose)
        return(final)

##### check for allele flips before we return results (when verbose == TRUE) #####
    if(nburn + niter > 0)
    {
        # how different are they from the original estimates from Pm.counts
        P.resid <- final$P - P.ori

        # don't count those markers that were unknown from the start -- they will have frequency of exactly 0.5
        P.resid <- ifelse(P.ori == 0.5, NA, P.resid)

        mu <- apply(P.resid, 2, mean, na.rm = TRUE)
        sigma <- apply(P.resid, 2, sd, na.rm = TRUE)

        P.deviation <- P.resid / matrix(rep(sigma, each = dim(P.resid)[1]), nrow = dim(P.resid)[1])

        ## if(any(abs(P.deviation) >= 2))
        ## {
        ##     bad.priors <- apply(abs(P.deviation) >= 2, 2, which)

        ##     flipped.deviation <- list()
        ##     if(length(bad.priors) > 0)
        ##     {
        ##         for(i in 1:length(bad.priors))
        ##             flipped.deviation[[names(bad.priors)[i]]] <- (final$P[bad.priors[[i]], i] -
        ##                                                          (1 - P.ori[bad.priors[[i]],i])) / sigma[i]

        ##         flips <- lapply(flipped.deviation, function(x) names(which(abs(x) < 2)))
        ##         bad.priors <- lapply(flipped.deviation, function(x) names(which(abs(x) >= 2)))

        ##         warning(length(unique(unlist(flips))), " potential flipped alleles and ",
        ##                 length(unique(unlist(bad.priors))), " additional outliers detected in Pm.counts.")
        ##     }
        ## }
    ## }else{
        bad.priors <- NULL
        flips <- NULL
    }

##### Include checks in retval (verbose == TRUE) #####

    ## final$sex.check <- sex.check
    ## final$indiv.check <- indiv.check
    ## final$marker.check <- marker.check
    ## final$hwe.check <- hwe
    ## final$flips <- flips
    ## final$bad.priors <- bad.priors

    class(final) <- 'admixture'

    return(final)
}

admixture.setup <- function(...)
{
    args <- list(...)
    if(is.null(args$burn))
        args$burn <- 100
    if(is.null(args$iter))
        args$iter <- 200
    if(is.null(args$every))
        args$every <- 5
    if(is.null(args$lambda))
        args$lambda <- 6
    if(is.null(args$tau))
        args$tau <- 300
    if(is.null(args$sex.chr))
        args$sex.chr <- 23
    if(is.null(args$indiv.geno.check))
        args$indiv.geno.check <- 0.98
    if(is.null(args$marker.geno.check))
        args$marker.geno.check <- 0.98
    if(is.null(args$hwe.thresh))
        args$hwe.thresh <- 1e-4
    if(is.null(args$male.het.X))
        args$male.het.X <- 1
    if(is.null(args$cores))
        args$cores <- detectCores()
    if(is.null(args$dev))
        args$dev <- FALSE
    if(is.null(args$verbose))
        args$verbose <- TRUE
    if(is.null(args$debug))
        args$debug <- FALSE
    if(is.null(args$fast))
        args$fast <- FALSE

    haps = args$haps
    geno = args$geno
    gender = args$gender
    chr = args$chr
    pos = args$pos
    Pm.prior = args$Pm.prior
    burn = args$burn
    iter = args$iter
    every = args$every
    indiv.id = args$indiv.id
    marker.id = args$marker.id
    pop.id = args$pop.id
    lambda = args$lambda
    tau = args$tau
    omega = args$omega
    rand.seed = args$rand.seed
    cores = args$cores
    cl = args$cl
    dev = args$dev
    verbose = args$verbose
    debug = args$debug
    sex.chr = args$sex.chr
    indiv.geno.check = args$indiv.geno.check
    marker.geno.check = args$marker.geno.check
    hwe.thresh = args$hwe.thresh
    bad.indiv = args$bad.indiv
    male.het.X = args$male.het.X
    fast = args$fast

    return(environment())
}
