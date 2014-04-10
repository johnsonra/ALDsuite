# admixture.R
# Infer ancestry estimates for a data set
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created September 10, 2012
# Last Modified December 11, 2013

# see documentation for details on these arguments
admixture <- function(Pm.prior, haps = NULL, geno = NULL, gender = NULL, chr, pos, burn = 100,
                      iter = 200, every = 5, indiv.id = NULL, marker.id = NULL, pop.id = NULL,
                      lambda = 9, tau = 300, omega = NULL, rand.seed = NULL,
                      cores = detectCores(), cl = NULL, dev = FALSE, verbose = TRUE,
                      debug = FALSE, sex.chr = 23, indiv.geno.check = 0.98,
                      marker.geno.check = .98, hwe.thresh = 1e-4, bad.indiv = NULL,
                      bad.marker = NULL, male.het.X = 1, fast = FALSE, A0 = NULL, Ak = NULL)
{
##########
# Checks #
##########

    #<<<<<<<<<<<<<<<<<< checks need to be updated to handle both haps and geno >>>>>>>>>>>>>>>>>

##### be sure all variables are formatted properly #####

    ## # type checks on ids
    ## if(is.factor(indiv.id))
    ##     indiv.id <- as.character(indiv.id)

    ## if(is.factor(marker.id))
    ##     marker.id <- as.character(marker.id)

    ## if(is.factor(pop.id))
    ##     pop.id <- as.character(pop.id)

    ## # boolean checks
    ## if(!is.logical(dev))
    ## {
    ##     warning("Incorrect data type for dev, ignoring input.")
    ##     dev <- FALSE
    ## }

    ## if(!is.logical(verbose))
    ## {
    ##     warning("Incorrect data type for verbose, ignoring input.")
    ##     verbose <- FALSE
    ## }

    ## if(!is.logical(debug))
    ## {
    ##     warning("Incorrect data type for debug, ignoring input.")
    ##     debug <- FALSE
    ## }

    ## # geno format check
    ## if(any(dimnames(geno)[[1]] != indiv.id))
    ##     warning(sum(dimnames(geno)[[1]] != indiv.id), " individual names in geno do not match indiv.id.")

    ## if(any(dimnames(geno)[[2]] != marker.id))
    ##     warning(sum(dimnames(geno)[[2]] != marker.id), " marker names in geno do not match marker.id.")

    ## # gender format check
    ## if(!is.null(gender) & length(gender) != dim(geno)[1])
    ##     stop("Length of gender should be ", dim(geno)[1], '.')

    ## if(any(names(gender) != indiv.id))
    ##     warning(sum(names(gender) != indiv.id), " names of the gender variable do not match indiv.id.")

    ## # chr format check
    ## if(length(chr) != length(Pm.prior))
    ##     stop("Length of chr should be ", length(Pm.prior), '.')

    ## if(any(names(chr) != marker.id))
    ##     warning(sum(names(chr) != marker.id), " names of the chr variable do not match marker.id.")

    ## # pos/d format check
    ## if(!is.null(pos) & length(pos) != dim(geno)[2])
    ##     stop("Length of pos should be", dim(geno)[2], '.')

    # Pm.prior format check
    ## if(!all(names(Pm.prior[[1]][[3]]) == pop.id) & !is.null(pop.id))
    ##     stop("Pm.prior must have the same names as pop.id when pop.id is not null.")

    ## names of Pm.prior should be in column names of haps (and geno?)
    ## if(!length(Pm.prior) == dim(geno)[2])
    ##     stop("Each element of Pm.prior must have the same length, equal to the number of columns in geno.")

    ## if(any(names(Pm.prior) != marker.id))
    ##     stop("Each marker in Pm.prior needs to be present in marker.id when not null.")

    ## if(any(unlist(lapply(Pm.prior, length)) != 3) | # each marker should have 3 elements
    ##    length(table(unlist(lapply(Pm.prior, lapply, length)))) != 1) # each element should have length K
    ##     stop("Formatting error in Pm.prior detected.")


    ## # numeric checks
    ## if(!is.numeric(burn) | length(burn) != 1 | burn < 0)
    ## {
    ##     warning("Invalid value for burn, ignoring input.")
    ##     burn <- 100
    ## }

    ## if(!is.numeric(iter) | length(iter) != 1 | iter < 0)
    ## {
    ##     warning("Invalid value for iter, ignoring input.")
    ##     iter <- 200
    ## }

    ## if(!is.numeric(every) | length(every) != 1 | every < 1)
    ## {
    ##     warning("Invalid value for every, ignoring input.")
    ##     every <- 5
    ## }

    ## if(!is.numeric(lambda) | length(lambda) != 1 | lambda <= 0)
    ## {
    ##     warning("Invalid value for lambda, ignoring input.")
    ##     lambda <- 6
    ## }

    ## if(!is.numeric(tau) | !length(tau) %in% c(1, length(pop.id))  | tau <= 0)
    ## {
    ##     warning("Invalid value for tau, ignoring input.")
    ##     tau <- 300
    ## }

    ## if(!is.null(omega) & (!is.numeric(omega) | any(omega < 0) |
    ##                       (length(omega) != length(pop.id) & !is.null(pop.id))))
    ## {
    ##     warning("Invalid value for omega, ignoring input.")
    ##     omega <- NULL
    ## }

    ## if(!is.null(rand.seed) & !is.numeric(rand.seed))
    ## {
    ##     warning("Invalid value for rand.seed, ignoring input.")
    ##     rand.seed <- NULL
    ## }else{
    ##     if(!is.null(rand.seed))
    ##         set.seed(rand.seed)
    ## }

    ## if(!is.numeric(cores) | cores < 1)
    ## {
    ##     warning("Invalid value for cores, ignoring input.")
    ##     cores <- detectCores()
    ## }

    ## if(!is.numeric(indiv.geno.check) | indiv.geno.check < 0 | indiv.geno.check > 1)
    ## {
    ##     warning("Bad value for indiv.geno.check, ignoring input.")
    ##     indiv.geno.check <- 0.98
    ## }

    ## if(!is.numeric(marker.geno.check) | marker.geno.check < 0 | marker.geno.check > 1)
    ## {
    ##     warning("Bad value for marker.geno.check, ignoring input.")
    ##     marker.geno.check <- 0.98
    ## }

    ## # IDs - some of this is already taken care of above
    ## if(!is.null(indiv.id) & length(indiv.id) != dim(geno)[1])
    ##     stop("Incorrect length of indiv.id.")

    ## # names of Pm.prior should be the same as marker.ids (and ordered properly)
    ## if(!is.null(marker.id) & length(marker.id) != length(Pm.prior))
    ##     stop("Incorrect length of marker.id.")

    ## if(!is.null(pop.id) & length(pop.id) != length(Pm.prior[[1]]$freq))
    ##     stop("Incorrect length of pop.id.")

    ## # other data types
    ## if(!is.null(cl) & !"cluster" %in% class(cl))
    ##     stop("Invalid value for cl.")


##### be sure gender is specified if anything in chr == sex.chr #####
    ## if(is.null(gender) & any(chr == sex.chr))
    ##     stop("Sex chromosomes detected but gender is null")


##### drop any bad individuals / makers #####

    #### make sure bad.marker and bad.indiv are reasonable before doing this!!!

    ## if(!is.null(bad.marker))
    ## {
    ##     warning("bad.marker is not properly tested. Odd things will happen when using this option!")
    ##     ## geno <- geno[,-bad.marker]
    ##     chr <- chr[-bad.marker]
    ##     marker.id <- marker.id[-bad.marker]
    ##     pos <- pos[-bad.marker]
    ##     ## Pm.counts <- Pm.counts[-bad.marker,,] !!!!!!!!!!!!!!!!!!! um...Pm.prior is going to be different..
    ## }

    ## if(!is.null(bad.indiv))
    ## {
    ##     ## geno <- geno[-bad.indiv,]
    ##     gender <- gender[-bad.indiv]
    ##     indiv.id <- indiv.id[-bad.indiv]
    ## }


##### gender check #####
    ## if(any(chr == sex.chr))
    ## {
    ##     # check males
    ##     ## if(any(na.omit(unlist(geno[gender == 'M', chr == sex.chr])) == 1))
    ##     ## {
    ##     ##     # get number of genotypes with more than one allele
    ##     ##     sex.check <- apply(geno[gender == 'M', chr == sex.chr] == 1, 1, sum, na.rm = TRUE)
    ##     ##     sex.check <- sex.check[sex.check > male.het.X]

    ##     ##     # make heterozygous loci missing
    ##     ##     geno[gender == 'M', chr == sex.chr] <- apply(geno[gender == 'M', chr == sex.chr], 1:2,
    ##     ##                                                  function(x) ifelse(x == 1, NA, x))

    ##     ##     # throw a warning about how many we have above the threshold
    ##     ##     if(length(sex.check) > 0)
    ##     ##         warning(length(sex.check), " males have between ", min(sex.check), " and ",
    ##     ##                 max(sex.check), " heterozygote X-chromosome genotypes.")
    ##     ## }else{
    ##     ##     sex.check <- NULL
    ##     ## }

    ##     # check females -- given the non-random selection of SNPs (based on frequencies in races) and
    ##     #                  the fact that some individuals might have very high admixture percentages
    ##     #                  for some races, I'm not sure how exactly to do this for females...
    ##     # PLINK uses the inbreeding coefficient: F ... might consider wether this would work here.
    ##     # If we compare the autosome F with the X chromosome F, we might be onto something...
    ##     # not convinced just yet
    ##     ## autosomes <- apply(geno[gender == 'F', chr != sex.chr] > 1, 1, sum) / sum(chr != sex.chr)
    ##     ## sex.check <- apply(geno[gender == 'F', chr == sex.chr] > 1, 1, sum) / sum(chr == sex.chr)
    ## }else{
        sex.check <- NULL
    ## }

##### check for complete genotyping (by individual / by marker) #####
    ## indiv.check <- apply(!is.na(haps), 1, sum) / dim(haps)[2]
    ## marker.check <- apply(!is.na(haps), 2, sum) / dim(haps)[1]

    ## if(any(indiv.check < indiv.geno.check))
    ##     warning(sum(indiv.check < indiv.geno.check), " individuals with less than ",
    ##             round(indiv.geno.check * 100),
    ##             "% complete genotyping detected.")

    ## if(any(marker.check < marker.geno.check))
    ##     warning(sum(marker.check < marker.geno.check), " markers with less than ",
    ##             round(marker.geno.check * 100),
    ##             "% complete genotyping detected.")

##### check for an excess of monomorphic markers #####
    ## # allele frequencies
    ## tmp <- sapply(names(Pm.prior), function(x) mean(haps[,x] / 2, na.rm = TRUE) / 2)

    ## if(any(tmp == 0))
    ## {
    ##     warning(sum(tmp == 0), " monomorphic markers present. A large number of monomorphic",
    ##             " markers could indicate some problems with the data.")
    ## }

    #################################################### need to check on this

##### check HWE #####
    hwe <- NULL
    ## hwe <- apply(geno[,chr != sex.chr], 2, table)
    ## if(is.list(hwe))
    ## {
    ##     tmp <- matrix(0, nrow = 3, ncol = length(hwe),
    ##                   dimnames = list(c('0', '1', '2'), colnames(geno)[chr != sex.chr]))

    ##     tmp['0',] <- sapply(hwe, `[`, '0')
    ##     tmp['1',] <- sapply(hwe, `[`, '1')
    ##     tmp['2',] <- sapply(hwe, `[`, '2')

    ##     tmp[is.na(tmp)] <- 0

    ##     hwe <- tmp
    ## }

    ## hwe <- hwexact(hwe['0',], hwe['1',], hwe['2',])
    ## names(hwe) <- marker.id[!chr == sex.chr]

    ## # only do females on this set...can't really say anything with males
    ## if(any(chr == sex.chr) & any(gender == 'F'))
    ## {
    ##     hweX <- apply(geno[gender == 'F', chr == sex.chr], 2, table)
    ##     if(is.list(hweX))
    ##     {
    ##         tmp <- matrix(0, nrow = 3, ncol = length(hweX),
    ##                       dimnames = list(c('0', '1', '2'), colnames(geno)[chr == sex.chr]))

    ##         tmp['0',] <- sapply(hweX, `[`, '0')
    ##         tmp['1',] <- sapply(hweX, `[`, '1')
    ##         tmp['2',] <- sapply(hweX, `[`, '2')

    ##         tmp[is.na(tmp)] <- 0

    ##         hweX <- tmp
    ##     }

    ##     hweX <- hwexact(hweX['0',], hweX['1',], hweX['2',])
    ##     names(hweX) <- marker.id[chr == sex.chr]

    ##     hwe <- c(hwe, hweX)
    ## }

    ## if(any(hwe < hwe.thresh))
    ##     warning(sum(hwe < hwe.thresh), " markers failed HWE test.")

##### check for allele flips later on in the process #####


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

######### haplotypes #########

    if(!is.null(haps))
    {
        tmp <- array(0, dim = c(dim(haps)[1] / 2, dim(haps)[2], 2),
                     dimnames = list(dimnames(haps)[[1]][1:(dim(haps)[1]/2) * 2 - 1],
                                     dimnames(haps)[[2]], c('Mother', 'Father')))

        tmp[,,1] <- haps[1:dim(tmp)[1] * 2 - 1,]
        tmp[,,2] <- haps[1:dim(tmp)[1] * 2,]
        haps <- tmp

        geno <- apply(haps, 1:2, sum)
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

        invisible(clusterEvalQ(cl, library(mald)))

        clusterExport(cl, list('hmm', 'haps', 'geno', 'P', 'lambda', 'lambdaX', 'd', 'chr', 'A0', 'Ak', 'AX',
                               'gender', 'sex.chr', 'dev', 'omega', 'omegaX', 'W', 'M', 'alpha', 'alphaX',
                               'Pm', 'Pm.counts', 'Pm.prior', 'every', 'tau', 'verbose', 'O',
                               'burnin', 'burn.lik', 'iter.lik', 'sigmaA', 'sigmaL', 'debug'),
                      envir = environment())
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
            print(supercyc)

            if(cores > 1)
            {
                clusterExport(cl, list('supercyc'), envir = environment())
                invisible(clusterEvalQ(cl, eval(hmm)))
            }else{

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

            if(cores > 1 & supercyc != nburn) # this is for the *next burn-in* sample (i.e. skip the last cycle)
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
                mixing <- supercyc / nburn # This is for the next cycle -- first time is essentially 1/0

                clusterExport(cl, list('Ak.r', 'AX.r', 'omega.r', 'omegaX.r', 'lambda.r',
                                       'lambdaX.r', 'alpha.r', 'beta.r', 'alphaX.r', 'betaX.r',
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

                                 beta <- mixing*beta + (1 - mixing)*beta.r
                                 betaX <- mixing*betaX + (1 - mixing)*betaX.r
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
                invisible(clusterApplyLB(cl, 1:iter, function(supercyc) eval(hmm, envir = globalenv())))
            }else{
                invisible(clusterEvalQ(cl, for(supercyc in 1:niter) eval(hmm)))
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
                final$beta <- final$beta + chains[[i]]$beta
                final$alphaX <- final$alphaX + chains[[i]]$alphaX
                final$tau <- final$tau + chains[[i]]$tau
                final$sigmaA <- final$sigmaA + chains[[i]]$sigmaA
                final$sigmaL <- final$sigmaL + chains[[i]]$sigmaL
            }

        }else{
            # running things serially here
            for(supercyc in 1:(niter / every))
            {
                print(supercyc)

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
            divisor <- iter / every
        }else{
            divisor <- cores * niter / every
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
    bad.marker = args$bad.marker
    male.het.X = args$male.het.X
    fast = args$fast

    return(environment())
}
