# pcr.prior.R
# Function to set up Principal Component Regression prior coefficients
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created August 26, 2013
# Last Modified December 10, 2013

# snps = vector of rsnumbers in our sample
# pops = vector of prior population ids (i.e. 'CEU', 'YRI', ...)
# thresh = percentage of variation we desire the PCs to explain
## setup.prior.dev <- function(snps, pops, anchors = NULL, thresh = 0.8, maxpcs = NULL, cM.linked = 0.1,
##                         window = 0.1, phased = FALSE, dev = FALSE, n.samp = 100)
## {
##     data(hapmap, envir = environment())

##     hapmap <- subset(hapmap, rs %in% snps)
##     hapmap <- hapmap[order(hapmap$chr, hapmap$pos),] # should already be ordered, but just in case...

##     # if anchors are not specified, choose them
##     if(is.null(anchors))
##     {
##         anchors <- character()

##         ## score <- function(x) apply(cbind(abs(x*log(x)), 0), 1, max, na.rm = TRUE)
##         score <- function(x) -(x - .5)^2 + .25

##         # calculate aims/selection score
##         hapmap$score <- 0

##         for(k1 in 1:(length(pops) - 1)) # comparisons with self aren't helpful
##         {
##             fa <- paste('f', tolower(pops[k1]), sep = '.')
##             for(k2 in (k1+1):length(pops))
##             {
##                 fb <- paste('f', tolower(pops[k2]), sep = '.')
##                 hapmap$score <- hapmap$score + abs(hapmap[[fa]] - hapmap[[fb]]) *
##                                ( score(hapmap[[fa]]) + score(hapmap[[fb]]) )
##             }
##         }


##         # drop any uninformative markers
##         hapmap <- subset(hapmap, score > 0)


##         # pack the map with spacing of about at least cM.linked cM
##         for(c in unique(hapmap$chr))
##         {
##             candidates <- subset(hapmap, chr == c)
##             while(dim(candidates)[1] > 0)
##             {
##                 pick <- candidates$rs[which.max(candidates$score)]

##                 candidates <- subset(candidates, abs(cM - subset(candidates, rs == pick)$cM) > window)

##                 anchors <- c(anchors, pick)
##             }
##         }

##         # reorder anchors correctly
##         anchors <- hapmap$rs[hapmap$rs %in% anchors]
##     }

##     # make general framework for Pm.prior
##     retval <- lapply(anchors, function(x) list(freq = numeric(),
##                                                n = numeric()))
##     names(retval) <- anchors


##     # generalize another time!
##     data(yri20, envir = environment())
##     yri <- phased
##     LD1 <- LD

##     data(ceu20, envir = environment())
##     ceu <- phased
##     LD2 <- LD

##     for(rs in anchors)
##     {
##         # setup freq and n
##         retval[[rs]]$freq['YRI'] <- mean(yri[,rs])
##         retval[[rs]]$freq['CEU'] <- mean(ceu[,rs])
##         retval[[rs]]$n['YRI'] <- dim(yri)[1]
##         retval[[rs]]$n['CEU'] <- dim(ceu)[1]

##         # pick linked alleles
##         hiLD <- unique(c(subset(LD1, rs1 == rs)$rs2, subset(LD1, rs2 == rs)$rs1,
##                          subset(LD2, rs1 == rs)$rs2, subset(LD2, rs2 == rs)$rs1))
##         loc <- hapmap$cM[hapmap$rs == rs]
##         linked <- unique(c(hiLD, subset(hapmap, abs(cM - loc) <= window & cM != loc)$rs))
##         ## linked <- subset(hapmap, abs(cM - loc) <= window & cM != loc & !rs %in% hiLD)$rs
##         linked <- linked[linked %in% hapmap$rs]

##         # set up model
##         tmp <- range(subset(hapmap, rs %in% c(linked) | cM == loc)$cM)
##         d <- tmp[2] - tmp[1]

##         retval[[rs]]$model <- list(YRI = list(linked = linked,
##                                               d = d,
##                                               eig = NULL,
##                                               betas = NULL,
##                                               vcv = NULL,
##                                               hessian = NULL),
##                                    CEU = list(linked = linked,
##                                               d = d,
##                                               eig = NULL,
##                                               betas = NULL,
##                                               vcv = NULL,
##                                               hessian = NULL))

##         # do pca of ceu and yri
##         if(length(linked) == 1)
##         {
##             lreg <- logitreg(yri[,rs], matrix(yri[,linked]))
##             retval[[rs]]$model[['YRI']]$eig <- 1
##             retval[[rs]]$model[['YRI']]$betas <- lreg$par
##             retval[[rs]]$model[['YRI']]$vcv <- lreg$vcv
##             retval[[rs]]$model[['YRI']]$hessian <- lreg$hessian

##             if(any(diag(retval[[rs]]$model$YRI$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$YRI$vcv)))
##             {
##                 retval[[rs]]$model$YRI$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$YRI$hessian)))
##                 retval[[rs]]$model$YRI$vcv <- solve(retval[[rs]]$model$YRI$hessian)
##             }

##             lreg <- logitreg(ceu[,rs], matrix(ceu[,linked]))
##             retval[[rs]]$model[['CEU']]$eig <- 1
##             retval[[rs]]$model[['CEU']]$betas <- lreg$par
##             retval[[rs]]$model[['CEU']]$vcv <- lreg$vcv
##             retval[[rs]]$model[['CEU']]$hessian <- lreg$hessian

##             if(any(diag(retval[[rs]]$model$CEU$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$CEU$vcv)))
##             {
##                 retval[[rs]]$model$CEU$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$CEU$hessian)))
##                 retval[[rs]]$model$CEU$vcv <- solve(retval[[rs]]$model$CEU$hessian)
##             }

##         }
##         if(length(linked) > 1)
##         {
##             # PCA
##             cormat <- cor(rbind(ceu[,linked], yri[,linked]), use = 'pairwise.complete.obs')
##             cormat[!is.finite(cormat)] <- 0
##             eig <- eigen(cormat)

##             if(sum(eig$values) == 0)
##             {
##                 retval[[rs]]$model$YRI$linked <- character()
##                 retval[[rs]]$model$CEU$linked <- character()
##             }else{
##                 # pick PCs
##                 npcs <- min(c(sum(cumsum(eig$values) / sum(eig$values) < thresh) + 1,
##                               length(eig$values), maxpcs))

##                 retval[[rs]]$model$YRI$eig <- eig$vectors[,1:npcs]
##                 retval[[rs]]$model$CEU$eig <- eig$vectors[,1:npcs]

##                 # Model PCs
##                 pcs <- yri[,linked] %*% retval[[rs]]$model$YRI$eig

##                 lreg <- logitreg(yri[,rs], pcs)
##                 retval[[rs]]$model$YRI$betas <- lreg$par
##                 retval[[rs]]$model$YRI$vcv <- lreg$vcv
##                 retval[[rs]]$model$YRI$hessian <- lreg$hessian

##                 if(any(diag(retval[[rs]]$model$YRI$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$YRI$vcv)))
##                 {
##                     retval[[rs]]$model$YRI$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$YRI$hessian)))
##                     retval[[rs]]$model$YRI$vcv <- solve(retval[[rs]]$model$YRI$hessian)
##                 }


##                 pcs <- ceu[,linked] %*% retval[[rs]]$model$CEU$eig

##                 lreg <- logitreg(ceu[,rs], pcs)
##                 retval[[rs]]$model$CEU$betas <- lreg$par
##                 retval[[rs]]$model$CEU$vcv <- lreg$vcv
##                 retval[[rs]]$model$CEU$hessian <- lreg$hessian

##                 if(any(diag(retval[[rs]]$model$CEU$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$CEU$vcv)))
##                 {
##                     retval[[rs]]$model$CEU$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$CEU$hessian)))
##                     retval[[rs]]$model$CEU$vcv <- solve(retval[[rs]]$model$CEU$hessian)
##                 }
##             }
##         }
##     }

##     return(retval)
## }

setup.prior.dev <- function(snps, pops, anchors = NULL, thresh = 0.8, maxpcs = NULL, cM.linked = 0.1,
                        window = 0.1, phased = FALSE, dev = FALSE, n.samp = 100)
{
    data(hapmap, envir = environment())

    hapmap <- subset(hapmap, rs %in% snps)
    hapmap <- hapmap[order(hapmap$chr, hapmap$pos),] # should already be ordered, but just in case...

    # if anchors are not specified, choose them
    if(is.null(anchors))
    {
        anchors <- character()

        ## score <- function(x) apply(cbind(abs(x*log(x)), 0), 1, max, na.rm = TRUE)
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
                               ( score(hapmap[[fa]]) + score(hapmap[[fb]]) )
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

    # make general framework for Pm.prior
    retval <- lapply(anchors, function(x) list(freq = numeric(),
                                               n = numeric()))
    names(retval) <- anchors


    # generalize another time!
    data(yri20, envir = environment())
    yri <- phased
    LD1 <- LD

    data(ceu20, envir = environment())
    ceu <- phased
    LD2 <- LD

    for(rs in anchors)
    {
        # setup freq and n
        retval[[rs]]$freq['YRI'] <- mean(yri[,rs])
        retval[[rs]]$freq['CEU'] <- mean(ceu[,rs])
        retval[[rs]]$n['YRI'] <- dim(yri)[1]
        retval[[rs]]$n['CEU'] <- dim(ceu)[1]

        # pick linked alleles
        hiLD <- unique(c(subset(LD1, rs1 == rs)$rs2, subset(LD1, rs2 == rs)$rs1,
                         subset(LD2, rs1 == rs)$rs2, subset(LD2, rs2 == rs)$rs1))
        loc <- hapmap$cM[hapmap$rs == rs]
        linked <- unique(c(hiLD, subset(hapmap, abs(cM - loc) <= window)$rs))
        ## linked <- subset(hapmap, abs(cM - loc) <= window & cM != loc & !rs %in% hiLD)$rs
        linked <- linked[linked %in% hapmap$rs]

        # set up model
        tmp <- range(subset(hapmap, rs %in% c(linked) | cM == loc)$cM)
        d <- tmp[2] - tmp[1]

        retval[[rs]]$model <- list(YRI = list(linked = linked,
                                              d = d,
                                              eig = NULL,
                                              betas = NULL,
                                              vcv = NULL,
                                              hessian = NULL),
                                   CEU = list(linked = character(),
                                              d = 0,
                                              eig = NULL,
                                              betas = NULL,
                                              vcv = NULL,
                                              hessian = NULL))

        # do pca of ceu and yri
        if(length(linked) == 1)
        {
            lreg <- logitreg(yri[,rs], matrix(yri[,linked]))
            retval[[rs]]$model[['YRI']]$eig <- 1
            retval[[rs]]$model[['YRI']]$betas <- lreg$par
            retval[[rs]]$model[['YRI']]$vcv <- lreg$vcv
            retval[[rs]]$model[['YRI']]$hessian <- lreg$hessian

            if(any(diag(retval[[rs]]$model$YRI$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$YRI$vcv)))
            {
                retval[[rs]]$model$YRI$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$YRI$hessian)))
                retval[[rs]]$model$YRI$vcv <- solve(retval[[rs]]$model$YRI$hessian)
            }

            lreg <- logitreg(ceu[,rs], matrix(ceu[,linked]))
            retval[[rs]]$model[['CEU']]$eig <- 1
            retval[[rs]]$model[['CEU']]$betas <- lreg$par
            retval[[rs]]$model[['CEU']]$vcv <- lreg$vcv
            retval[[rs]]$model[['CEU']]$hessian <- lreg$hessian

            if(any(diag(retval[[rs]]$model$CEU$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$CEU$vcv)))
            {
                retval[[rs]]$model$CEU$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$CEU$hessian)))
                retval[[rs]]$model$CEU$vcv <- solve(retval[[rs]]$model$CEU$hessian)
            }

        }else{
            # PCA
            cormat <- cor(rbind(ceu[,linked], yri[,linked]), use = 'pairwise.complete.obs')
            cormat[!is.finite(cormat)] <- 0
            eig <- eigen(cormat)

            if(sum(eig$values) == 0)
            {
                retval[[rs]]$model$YRI$linked <- character()
            }else{
                # pick PCs
                npcs <- min(c(sum(cumsum(eig$values) / sum(eig$values) < thresh) + 1,
                              length(eig$values), maxpcs))

                retval[[rs]]$model$YRI$eig <- eig$vectors[,1:npcs]

                # Model PCs
                yri.pcs <- yri[,linked] %*% retval[[rs]]$model$YRI$eig
                ceu.pcs <- ceu[,linked] %*% retval[[rs]]$model$YRI$eig

                y <- c(rep(1, dim(yri)[1]), rep(0, dim(ceu)[1]))
                x <- rbind(yri.pcs, ceu.pcs)

                lreg <- logitreg(y, x)
                retval[[rs]]$model$YRI$betas <- lreg$par
                retval[[rs]]$model$YRI$vcv <- lreg$vcv
                retval[[rs]]$model$YRI$hessian <- lreg$hessian

                if(any(diag(retval[[rs]]$model$YRI$vcv) <= 0) | any(!is.finite(retval[[rs]]$model$YRI$vcv)))
                {
                    retval[[rs]]$model$YRI$hessian <- rwish(2, diag(nrow = nrow(retval[[rs]]$model$YRI$hessian)))
                    retval[[rs]]$model$YRI$vcv <- solve(retval[[rs]]$model$YRI$hessian)
                }
            }
        }
    }

    return(retval)
}

setup.prior <- function(snps, pops, anchors = NULL, thresh = 0.8, maxpcs = 6, cM.linked = 0.1,
                        window = 0.1, phased = FALSE, dev = FALSE, n.samp = 100)
{
    data(hapmap, envir = environment())

    hapmap <- subset(hapmap, rs %in% snps)
    hapmap <- hapmap[order(hapmap$chr, hapmap$pos),] # should already be ordered, but just in case...

    # if anchors are not specified, choose them
    if(is.null(anchors))
    {
        anchors <- character()

        ## score <- function(x) apply(cbind(abs(x*log(x)), 0), 1, max, na.rm = TRUE)
        score <- function(x) -(x - .5)^2 + .25

        # calculate aims/selection score
        hapmap$score <- 0

        if(dev)
        {
            print('picking common markers')
            for(k in 1:length(pops))
            {
                f <- hapmap[[paste('f', tolower(pops[[k]]), sep = '.')]]
                f <- f + ifelse(f > 0.5, -0.001, 0.001)

                hapmap$score <- hapmap$score + log(f) + log(1 - f)
            }
        }else{
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
        }

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

    if(phased)
    {
        for(k in pops)
        {
            for(i in unique(hapmap$chr))
            {
                hapmap.sub <- subset(hapmap, chr == i)
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

                for(rs in hapmap.sub$rs[hapmap.sub$rs %in% anchors])
                    retval[[rs]] <- pcr.prior.phased(retval[[rs]], rs, toupper(k), hapmap.sub, phased, LD,
                                                     thresh, maxpcs, cM.linked, window, dev)

            }
        }

        class(retval) <- 'phased'

    }else{
        for(k1 in 1:length(pops))
        {
            for(k2 in k1:length(pops))
            {
                print(paste(k1, k2))
                for(c in unique(hapmap$chr))
                {
                    hapmap.sub <- subset(hapmap, chr == c)

                    # load
                    tmp1 <- load.pop(c, pops[k1], hapmap.sub$rs)
                    if(k1 != k2)
                    {
                        tmp2 <- load.pop(c, pops[k2], hapmap.sub$rs)
                    }else{
                        tmp2 <- NULL
                    }

                    for(rs in hapmap.sub$rs[hapmap.sub$rs %in% anchors])
                        retval[[rs]] <- pcr.prior.unphased(retval[[rs]], rs, toupper(c(pops[k1], pops[k2])),
                                                  hapmap.sub, list(tmp1, tmp2), thresh, maxpcs, cM.linked,
                                                  window, dev, n.samp)
                }
            }
        }
        class(retval) <- 'unphased'
    }

    class(retval) <- c('Pm.prior', class(retval))

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

# rs = the SNP id we are calculating the priors for
# map = subset of hapmap containing all sampled SNPs
# phased = correct set of phased chromosomes
# LD = correct set of LD measures
# thresh = percent of variation we want the PCs to explain
pcr.prior.phased <- function(retval, rs.cur, pop, map, phased, LD, thresh, maxpcs, cM.linked, window, dev)
{
    retval$freq[pop] <- mean(phased[,rs.cur])
    retval$n[pop] <- dim(phased)[1]

    # figure out what is linked
    hiLD <- character()#c(subset(LD, rs2 == rs.cur)$rs1, subset(LD, rs1 == rs.cur)$rs2)

    linked <- unique(c(subset(LD, rs2 == rs.cur)$rs1,
                       subset(LD, rs1 == rs.cur)$rs2,
                       subset(map, abs(cM - map$cM[map$rs == rs.cur]) <= window &
                                   !rs %in% c(rs.cur, hiLD))$rs))


    # drop markers outside of cM.linked
    ## drop <- subset(map, abs(cM - map$cM[map$rs == rs.cur]) > cM.linked & rs %in% linked)
    ## linked <- linked[!linked %in% drop]

    tmp <- range(subset(map, rs %in% c(rs.cur, linked))$cM)

    # model
    model <- list(linked = linked,
                  d = tmp[2] - tmp[1], # size of block in cM
                  eig = NULL,
                  betas = NULL,
                  vcv = NULL,
                  hessian = NULL)

    if(length(linked) == 0)
    {
        retval$model[[pop]] <- model
        return(retval)
    }

    # if only one marker, no need to calculate eigenvectors
    if(length(linked) == 1)
    {
        lreg <- logitreg(phased[,rs.cur], matrix(phased[,linked]))

        model$eig <- 1
        model$betas <- lreg$par
        model$vcv <- lreg$vcv
        model$hessian <- lreg$hessian
    }else{
        # calculate eigenvectors
        phased.sub <- phased[,linked]

        cormat <- cor(phased.sub, use='pairwise.complete.obs')
        cormat[!is.finite(cormat)] <- 0
        eig <- eigen(cormat)

        # This happens sometimes when all the "linked" markers are monomorphic in one population
        # These are presumably common in other populations, though, so keep them!
        if(sum(eig$values) == 0)
        {
            model$linked <- character()
            retval$model[[pop]] <- model
            return(retval)
        }

        # pick PCs
        npcs <- min(c(sum(cumsum(eig$values) / sum(eig$values) < thresh) + 1, maxpcs))

        model$eig <- eig$vectors[,1:npcs]

        # calculate principle components
        pcs <- phased.sub %*% model$eig

        # calculate betas for this marker
        lreg <- logitreg(phased[,rs.cur], pcs)

        # record results
        model$betas <- lreg$par
        model$vcv <- lreg$vcv
        model$hessian <- lreg$hessian
    }

    # some of these can be messed up...if any of the diagonals end up being nonpositive,
    # choose random sample - variance must be positive!!!
    if(any(diag(model$vcv) <= 0) | any(!is.finite(model$vcv)))
    {
        model$hessian <- rwish(2, diag(nrow = nrow(model$hessian)))
        model$vcv <- solve(model$hessian)
    }

    retval$model[[pop]] <- model

    return(retval)
}

pcr.prior.unphased <- function(retval, rs.cur, pops, map, dat, thresh, maxpcs, cM.linked, window, dev, n.samp)
{
    cM.linked <- max(cM.linked, window)

    # all pops will have a chance to be calculated here
    retval$freq[pops[1]] <- mean(dat[[1]]$phased[,rs.cur])
    retval$n[pops[1]] <- dim(dat[[1]]$phased)[1]

    # figure out what is linked
    linked <- with(dat[[1]], c(subset(LD, rs2 == rs.cur)$rs1,
                               subset(LD, rs1 == rs.cur)$rs2,
                               subset(map, abs(cM - map$cM[map$rs == rs.cur]) <= cM.linked &
                                      rs != rs.cur)$rs))
    if(!is.null(dat[[2]])) # might be null because we could have YRI/YRI for example
        linked <- c(linked, with(dat[[2]], c(subset(LD, rs2 == rs.cur)$rs1,
                                             subset(LD, rs1 == rs.cur)$rs2,
                                             subset(map, abs(cM - map$cM[map$rs == rs.cur]) <= cM.linked &
                                                    rs != rs.cur)$rs)))
    linked <- unique(linked)

    tmp <- range(subset(map, rs %in% c(rs.cur, linked))$cM)

    # model
    model <- list(linked = linked,
                  d = tmp[2] - tmp[1], # size of block in cM
                  eig = NULL,
                  betas = NULL,
                  vcv = NULL,
                  hessian = NULL)

    if(length(linked) == 0)
    {
        retval$model[[paste(pops[1], pops[2], sep = '.')]] <- model
        return(retval)
    }

    # sample genomes
    genos <- matrix(nrow = n.samp, ncol = length(linked) + 1)
    na <- dim(dat[[1]]$phased)[1]
    nb <- dim(dat[[2]]$phased)[1]

    for(i in 1:n.samp)
    {
        a <- dat[[1]]$phased[sample(1:na, size = 1), c(rs.cur, linked)]
        if(is.null(nb)) # homozygous for popA
        {
            b <- dat[[1]]$phased[sample(1:na, size = 1), c(rs.cur, linked)]
        }else{
            b <- dat[[2]]$phased[sample(1:nb, size = 1), c(rs.cur, linked)]
        }

        genos[i,] <- a + b
    }

    # if only one marker, no need to calculate eigenvectors
    if(length(linked) == 1)
    {
        lreg <- multilogitreg(cbind(genos[,1] == 0, genos[,1] == 1, genos[,1] == 2),
                              cbind(as.numeric(genos[,2] == 1), as.numeric(genos[,2] == 2)))

        model$eig <- diag(2)
        model$betas <- lreg$par
        model$vcv <- lreg$vcv
        model$hessian <- lreg$hessian
    }else{
        # calculate eigenvectors
        genos.sub <- cbind(genos[,-1] == 1, genos[,-1] == 2)

        cormat <- cor(genos.sub, use='pairwise.complete.obs')
        cormat[!is.finite(cormat)] <- 0
        eig <- eigen(cormat)

        # This happens sometimes when all the "linked" markers are monomorphic in one population
        # These are presumably common in other populations, though, so keep them!
        if(sum(eig$values) == 0)
        {
            model$linked <- character()
            retval$model[[paste(pops[1], pops[2], sep = '.')]] <- model
            return(retval)
        }

        # pick PCs
        npcs <- min(c(sum(cumsum(eig$values) / sum(eig$values) < thresh) + 1, maxpcs))

        model$eig <- eig$vectors[,1:npcs]

        # calculate principle components
        pcs <- genos.sub %*% model$eig

        # calculate betas for this marker
        lreg <- multilogitreg(cbind(genos[,1] == 0, genos[,1] == 1, genos[,1] == 2), pcs)

        # record results
        model$betas <- lreg$par
        model$vcv <- lreg$vcv
        model$hessian <- lreg$hessian
    }

    # some of these can be messed up...if any of the diagonals end up being nonpositive,
    # choose random sample - variance must be positive!!!
    if(any(diag(model$vcv) <= 0) | any(!is.finite(model$vcv)))
    {
        model$hessian <- rwish(2, diag(nrow = nrow(model$hessian)))
        model$vcv <- solve(model$hessian)
    }

    retval$model[[paste(pops[1], pops[2], sep = '.')]] <- model

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
                      dimnames = list(paste('b', 1:ncol(X) - 1, sep = ''), c('p1', 'p2')))

    if(hessian)
        fit$vcv <- solve.approx(fit$hessian) # variance covariance matrix

    invisible(fit)
}
