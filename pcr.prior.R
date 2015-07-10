# pcr.prior.R
# Function to set up Principal Component Regression prior coefficients
# Randall Johnson
# CCR Collaborative Bioinformatics Resource
# Leidos Biomedical Research, Inc

library(parallel)
library(foreach)
library(doMC)
registerDoMC(2)

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

setup.prior <- function(snps, pops, anchors = NULL, thresh = 0.8, maxpcs = 6,
                        window = 0.1, unphased = TRUE, n.samp = 300)
{
    if(!require(ALDdata))
        stop('ALDdata package required for setup.prior()')

    data(hapmap, envir = environment())

    hapmap <- subset(hapmap, rs %in% snps)
    hapmap <- hapmap[order(hapmap$chr, hapmap$pos),] # should already be ordered, but just in case...

    # if anchors are not specified, choose them
    if(is.null(anchors))
    {
        anchors <- character()

        score <- function(x) -(x - .5)^2 + .25

        # calculate aims/selection score
        hapmap$score <- 0

        for(k1 in 1:(length(pops) - 1)) # comparisons with self aren't helpful
        {
            fa <- paste('f', tolower(pops[k1]), sep = '.')
            for(k2 in (k1+1):length(pops))
            {
                fb <- paste('f', tolower(pops[k2]), sep = '.')
                hapmap$score <- hapmap$score + abs(hapmap[[fa]] - hapmap[[fb]]) *
                    (score(hapmap[[fa]]) + score(hapmap[[fb]]))
            }
        }


        # drop any uninformative markers
        hapmap <- subset(hapmap, score > 0)


        # pack the map with spacing of about at least cM.linked cM
        for(c in unique(hapmap$chr))
        {
            candidates <- subset(hapmap, chr == c)
            while(dim(candidates)[1] > 0)
            {
                pick <- candidates$rs[which.max(candidates$score)]

                candidates <- subset(candidates, abs(cM - subset(candidates, rs == pick)$cM) > window)

                anchors <- c(anchors, pick)
            }
        }

        # reorder anchors correctly
        anchors <- hapmap$rs[hapmap$rs %in% anchors]
    }

    # initiate retval with "proper" ordering
    retval <- lapply(anchors, function(x) list(freq = numeric(),
                                               n = numeric()))
    names(retval) <- anchors


    ######### Set up priors #########
    # this is where the phased haplotype training data will reside
    train <- list(dat = NULL)

    for(i in unique(hapmap$chr))
    {
        hapmap.sub <- subset(hapmap, chr == i)

        ### Collect Traning Data for Chromosome i ###
        for(k in pops)
        {
            if(k == 'CHB+JPT')
            {
                eval(parse(text = paste('data(chb', i, ', envir = environment())', sep = '')))
                tmp <- phased
                tmp2 <- LD

                eval(parse(text = paste('data(jpt', i, ', envir = environment())', sep = '')))
                phased <- rbind(tmp, phased)

                LD <- merge(tmp2, LD, all = TRUE)
            }else{
                eval(parse(text = paste('data(', tolower(k), i, ', envir = environment())', sep = '')))
            }

            phased <- phased[,hapmap.sub$rs]
            LD <- subset(LD, rs1 %in% colnames(phased) & rs2 %in% colnames(phased))

            train$dat <- rbind(train$dat, cbind(which(pops == k), phased[,hapmap.sub$rs]))
            train[[k]] <- subset(LD, rs1 %in% colnames(phased) & rs2 %in% colnames(phased))
        }


        foreach(rs %doPar% hapmap.sub$rs[hapmap.sub$rs %in% anchors])
        {
            ### Calculate PC Loadings ###
            tmp <- pca.setup(train, rs, hapmap.sub, thresh, maxpcs, window, n.samp, unphased)

            ### PC Regression ###
            retval[[rs]] <- pcr.prior(tmp)

            ### Drop training data ###
            retval[[rs]]$train <- NULL
        }
    }

    if(unphased)
    {
        class(retval) <- c('Pm.prior', 'unphased')
    }else{
        class(retval) <- c('Pm.prior', 'phased')
    }

    return(retval)
}

load.pop <- function(chr, k, rsList)
{
    if(k == 'CHB+JPT')
    {
        eval(parse(text = paste('data(chb', chr, ', envir = environment())', sep = '')))
        tmp <- phased
        tmp2 <- LD

        eval(parse(text = paste('data(jpt', chr, ', envir = environment())', sep = '')))
        phased <- rbind(tmp, phased)

        LD <- merge(tmp2, LD, all = TRUE)
    }else{
        eval(parse(text = paste('data(', tolower(k), chr, ', envir = environment())', sep = '')))
    }

    phased <- phased[,rsList]
    LD <- subset(LD, rs1 %in% colnames(phased) & rs2 %in% colnames(phased))

    return(list(phased = phased, LD = LD))
}

# train = training data (list with one element per population, each with phased and LD data)
# rs.cur = the SNP id we are calculating the priors for
# map = subset of hapmap containing all sampled SNPs
# thresh = percent of variation we want the PCs to explain
# maxpcs = maximum number of PCs to track
# window = window size (in cM)
pca.setup <- function(train, rs.cur, map, thresh, maxpcs, window, n.samp, unphased)
{
    pops <- names(train)[-1]

    # figure out what is linked in any population (including rs.cur)
    linked <- subset(map, abs(cM - map$cM[map$rs == rs.cur]) <= window)$rs # close by

    # add in any markers with high LD
    for(k in pops)
    {
        linked <- unique(c(linked, subset(train[[k]], rs2 == rs.cur)$rs1,
                                   subset(train[[k]], rs1 == rs.cur)$rs2))
    }

    linked <- linked[linked %in% map$rs] # just to be sure...

    # set up model list
    model <- list(freq = numeric(),
                  n = numeric(),
                  linked = linked,
                  d = diff(range(subset(map, rs %in% linked)$cM)), # size of block in cM
                  eig = NULL,
                  betas = NULL,
                  vcv = NULL,
                  hessian = NULL,
                  train = cbind(train$dat[,1], train$dat[,linked]))

    for(k in 1:length(pops))
    {
        model$freq[pops[k]] <- mean(train$dat[train$dat[,1] == k, rs.cur]) / 2
        model$n[pops[k]] <- sum(train$dat[,1] == k)
    }

    if(unphased)
    {
        # combine chromosomes to simulate new people (these will have no crossovers --
        #   i.e. one chromosome from each ancestral population -- or two from one population)
        genos <- matrix(nrow = 0, ncol = dim(model$train)[2],
                        dimnames = list(character(), colnames(model$train)))

        for(k1 in 1:length(pops))
        {
            for(k2 in k1:length(pops))
            {
                pick1 <- sample((1:dim(model$train)[1])[model$train[,1] == k1],
                                size = n.samp, replace = TRUE)
                pick2 <- sample((1:dim(model$train)[1])[model$train[,1] == k2],
                                size = n.samp, replace = TRUE)

                # calculate genotypes
                pick <- cbind(as.numeric(paste(k1, k2, sep = '.')),
                              train$dat[pick1,linked] + train$dat[pick2,linked])

                # add newly simulated individuals to the unphased data set
                genos <- rbind(genos, pick)
            }
        }

        # replace phased data with unphased data
        model$train <- genos
    }

    # if only one marker, no need to calculate eigenvectors
    if(length(linked) == 1)
    {
        model$eig <- 1
    }else{
        # PCA for all linked markers
        cormat <- cor(model$train[,-1], use = 'pairwise.complete.obs')
        cormat[!is.finite(cormat)] <- 0

        eig <- eigen(cormat)

        if(sum(eig$values) == 0)
        {
            stop('deal with this')
        }else{
            # pick PCs
            npcs <- min(c(sum(cumsum(eig$values) / sum(eig$values) < thresh) + 1,
                          length(eig$values), maxpcs))

            # save PCs
            model$eig <- eig$vectors[,1:npcs]
        }
    }

    return(model)
}

pcr.prior <- function(retval)
{
    # responses
    groups <- unique(retval$train[,1])
    pops <- names(retval$freq)

    # regression
    if(dim(retval$train)[2] == 2)
    {
        X <- as.matrix(retval$train[,-1])
    }else{
        X <- retval$train[,-1]
    }

    if(length(groups) == 2)
    {
        y <- retval$train[,1] == groups[1]

        lreg <- logitreg(y, X %*% retval$eig)

        retval$betas <- t(t(lreg$par)) # collect and name betas
        colnames(retval$betas) <- pops[1]
    }else{
        y <- sapply(groups, function(i) retval$train[,1] == i)

        lreg <- multilogitreg(y, X %*% retval$eig)

        retval$betas <- lreg$par # collect and name betas

        groups <- strsplit(as.character(groups), '.', fixed = TRUE)
        groups[[length(groups)]] <- NULL
        colnames(retval$betas) <- paste(pops[as.numeric(sapply(groups, `[`, 1))],
                                        pops[as.numeric(sapply(groups, `[`, 2))], sep = '.')
    }

    # collect variance/covariance and hessian matrices
    retval$vcv <- lreg$vcv
    retval$hessian <- lreg$hessian


    # some of these can be messed up...if any of the diagonals end up being nonpositive,
    # randomly sample - variance must be positive!!!
    if(any(diag(retval$vcv) <= 0) | any(!is.finite(retval$vcv)))
    {
        retval$hessian <- rwish(2, diag(nrow = nrow(retval$hessian)))
        retval$vcv <- solve(retval$hessian)
    }

    return(retval)
}

# I need this because the Fisher scoring algorithm used by default in glm
# has problems when there is a perfect fit of the data...this happens in
# our case when tightly linked markers are included! Is it a problem???
# see MASS book pages 197-8 and 445
logitreg <- function(y, x, wt = rep(1, length(y)), intercept = TRUE, start = rep(0, p), hessian = TRUE, ...)
{
    # if x has no columns, reutrn odds based only on y
    if(ncol(x) == 0)
        return(mean(y) / (1 - mean(y)))

    fmin <- function(beta, X, y, w)
    {
        p <- plogis(X %*% beta)
        -sum(2 * w * ifelse(y, log(p), log(1-p)))
    }

    gmin <- function(beta, X, y, w)
    {
        eta <- X %*% beta
        p <- plogis(eta)

        -2 * t(w * dlogis(eta) * ifelse(y, 1/p, -1 / (1 - p))) %*% X
    }

    if(is.null(dim(x)))
       dim(x) <- c(length(x), 1)

    dn <- dimnames(x)[[2]]

    if(!length(dn))
        dn <- paste("PC", 1:ncol(x), sep = "")

    p <- ncol(x) + intercept

    if(intercept)
    {
        x <- cbind(1, x)
        dn <- c("(Intercept)", dn)
    }

    if(is.factor(y))
        y <- unclass(y) != 1

    fit <- optim(start, fmin, gmin, X=x, y=y, w=wt, method = "BFGS", hessian = hessian, ...)

    names(fit$par) <- dn

    if(hessian)
        fit$vcv <- solve.approx(fit$hessian) # variance covariance matrix

    invisible(fit)
}

# adaptation of logitreg for multinomial logistic regression
# took out a lot from the function above, but used the relevant parts...
# see also http://en.wikipedia.org/wiki/Multinomial_logistic_regression
multilogitreg <- function(y, X, start = NULL, hessian = TRUE)
{
    if(!is.logical(y))
        stop('y must be a logical matrix')

    # add intercept
    X <- cbind(1, X)

    # if start is unspecified, make 0
    start <- rep(0, ncol(X) * (ncol(y) - 1))

    # multinomial logistic regression function to minimize
    fmin <- function(Beta, X, y)
    {
        # ratios
        lr <- X %*% matrix(Beta, ncol(X), ncol(y) - 1)

        # probability of last category
        lp0 <- -log(1 + apply(exp(lr), 1, sum))

        if(!all(lp0 > -Inf))
        {
            # summing with respect to the bigger part
            # see http://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation.2Fsubtraction
            # if this results in errors at some point, look at how I did something similar in P.O2
            bigger <- apply(lr, 1, max)
            smaller <- apply(lr, 1, min)

            lp0 <- -(bigger + log(1 + exp(log(1 + exp(smaller)) - bigger)))

            lp <- cbind(lp0, lp0 + lr)

        }else{
            # log probabilities
            lp <- cbind(lp0, lp0 + lr)
        }

        return(-sum(lp[as.logical(y)]))
    }

    # calculate the gradiant and this will go a lot faster!
    gmin <- function(Beta, X, y)
    {
        r <- exp(X %*% Beta)

        sum.r <- apply(r, 1, sum)

        lplast <- -log(1 + sum.r)

        dlplast <- -1 / (1 + sum.r)^2 * (X %*% r)
        dlp <- 1 / lplast


    }

    fit <- optim(start, fmin, X = X, y = y,# method = 'BFGS',
                 hessian = hessian)

    fit$par <- matrix(fit$par, ncol(X), ncol(y) - 1,
                      dimnames = list(c('(Intercept)', paste('PC', 1:(ncol(X) - 1), sep = '')), NULL))

    if(hessian)
        fit$vcv <- solve.approx(fit$hessian) # variance covariance matrix

    invisible(fit)
}