# scripts.R
# nice short snipits of code that will be useful

##### Set up environment for sample.gammas.marg #####
tmp <- new.env()
tmp$fun <- 'sample.G'
tmp$geno <- haps
tmp$P <- P
tmp$Pm.prior <- Pm.prior
tmp$lambda <- lambda
tmp$lambdaX <- lambdaX
tmp$d <- d
tmp$chr <- chr
tmp$A0 <- A0
tmp$Ak <- Ak
tmp$AX <- AX
tmp$gender <- gender
tmp$sex.chr <- sex.chr
tmp$haps <- TRUE
tmp$dev <- dev

tmp <- new.env()
tmp$fun <- 'sample.G'
tmp$geno <- geno
tmp$P <- P
tmp$Pm.prior <- Pm.prior.geno
tmp$lambda <- lambda
tmp$lambdaX <- lambdaX
tmp$d <- d
tmp$chr <- chr
tmp$A0 <- A0
tmp$Ak <- Ak
tmp$AX <- AX
tmp$gender <- gender
tmp$sex.chr <- sex.chr
tmp$haps <- FALSE
tmp$dev <- dev

##### logP.G.lambda.d #####
tmp <- new.env()
tmp$fun <- 'logP.G.lambda.d'
tmp$G <- G.new
tmp$lambda <- lambda
tmp$d <- d
tmp$Ak <- Ak

##### P.a.gamma #####
tmp <- new.env()
tmp$fun <- 'P.a.gamma'
tmp$locus <- geno[,j,]
tmp$geno <- geno
tmp$prior <- Pm.prior[[j]]$model
tmp$P <- P[j,]
tmp$Ak <- Ak
tmp$lambda <- lambda
tmp$d <- d[j]

tmp$prev <- gammas[,max(j-1, 1),,]

##### P.gammas.a #####
tmp <- new.env()
tmp$fun <- 'P.gammas.a'
tmp$locus <- geno[,j]
tmp$geno <- geno
tmp$prior <- Pm.prior[[j]]
tmp$P <- P[j,]
tmp$A0 <- A0
tmp$Ak <- Ak
tmp$lambda <- lambda
tmp$d <- d[j]
tmp$gammas.prev <- gammas[,max(j - 1, 1),1,]

##### P.O2 #####
tmp <- new.env()
tmp$fun <- 'P.O2'
tmp$O.prev <- geno[,prior$model[[k]]$linked, drop = FALSE]
tmp$P <- P[c(k1, k2)]
tmp$prior <- prior$model[[k]]
tmp$linked <- length(prior$model[[k]]$linked) > 0

##### sample.gammas #####
tmp <- new.env()
tmp$fun <-  'sample.gammas'
tmp$pag <- pag
tmp$lambda <- lambda
tmp$Ak <- Ak
tmp$chr <- chr
tmp$pos <- pos
tmp$supercyc <- supercyc
tmp$burn <- 50
tmp$iter <- 100
tmp$every <- 5
tmp$G <- G
tmp$r <- r

##### Set up environemnt for updt.P.betas #####
tmp <- new.env()
tmp$fun <- 'updt.P.betas'
tmp$P <- P
tmp$Pm <- Pm
tmp$Pm.prior <- Pm.prior
tmp$lambda <- lambda
tmp$tau <- tau
tmp$G <- G
if(is.null(haps)){
    tmp$haps <- geno
}else{
    tmp$haps <- haps}
tmp$phased <- !is.null(haps)
tmp$dev <- dev

##### Set up environment for updt.betas #####
tmp <- new.env()
tmp$fun <- 'updt.betas'
tmp$j <- j
tmp$snp <- names(Pm.prior)[j]
tmp$prior <- Pm.prior[[j]]$model[[k]]
tmp$G <- G
tmp$haps <- haps
tmp$k <- k
tmp$tau <- tau
tmp$phased <- phased
tmp$dev <- dev

##### Set up environment for setup.prior #####
tmp <- new.env()
tmp$fun <- 'setup.prior'
tmp$snps <- hapmap$rs
tmp$pops <- c('YRI', 'CEU')
tmp$thresh <- 0.8
tmp$maxpcs <- NULL
tmp$cM.linked <- 0.1
tmp$anchors <- NULL
tmp$dev <- TRUE
tmp$n.samp <- 100
tmp$window <- 0.1
tmp$pahsed <- TRUE

##### Set up environment for pcr.prior #####
tmp <- new.env()
tmp$fun <- 'pcr.prior'
tmp$retval <- retval[[rs]]
tmp$rs.cur <- rs
tmp$pops <- toupper(c(pops[k1], pops[k2]))
tmp$map <- hapmap.sub
tmp$dat <- list(tmp1, tmp2)
tmp$thresh <- thresh
tmp$macpcs <- maxpcs
tmp$cM.linked <- cM.linked
tmp$window <- window
tmp$dev <- dev
tmp$n.samp <- n.samp

##### Set up environment for pcr.prior.one.way #####
tmp3 <- new.env()
tmp3$rs <- rs
tmp3$linked <- subset(LD, rs2 == rs)$rs1
tmp3$phased <- phased
tmp3$thresh <- thresh

##### Set up environment for logitreg #####
tmp4 <- new.env()
tmp4$y <- phased[,rs]
tmp4$X <- pcs
tmp4$wt <- rep(1, length(tmp4$y))
tmp4$intercept <- TRUE
tmp4$p <- ncol(as.matrix(phased[,rs])) + 1
tmp4$start <- rep(0, tmp4$p)

##### Set up environment for multilogitreg #####
tmp <- new.env()
tmp$fun <- 'multilogitreg'
tmp$y <- cbind(genos[,1] == 0, genos[,1] == 1, genos[,1] == 2)
tmp$X <- pcs
tmp$start <- NULL
tmp$hessian <- TRUE

##### Set up environment for sample.G #####
tmp2 <- new.env()
tmp2$pag <- pag
tmp2$lambda <- lambda
tmp2$Ak <- Ak
tmp2$chr <- chr
tmp2$pos <- pos
tmp2$burn <- 50
tmp2$iter <- 100
tmp2$every <- 5

## # old
## tmp2 <- new.env()
## tmp2$gammas.joint <- gammas.joint
## tmp2$O <- O
## tmp2$pops <- dimnames(P)[[2]]
## tmp2$chr <- chr
## tmp2$gender <- gender
## tmp2$sex.chr <- sex.chr

##### Set up environment for sample.O #####
tmp2 <- new.env()
tmp2$a <- geno[,j]
tmp2$pO <- pO
tmp2$gammas <- ps

##### Set up environment for P.O #####
tmp2 <- new.env()
tmp2$fun <- 'P.O'
tmp2$O.prev <- haps[,Pm.prior[[j]]$model[[k]]$linked, drop = FALSE]
tmp2$prior <- Pm.prior[[j]]$model[[k]]
tmp2$linked <- length(Pm.prior[[j]]$model[[k]]$linked) > 0
tmp2$P <- P[j,][c(as.numeric(substr(combos[k], 2, 2)), as.numeric(substr(combos[k], 3, 3)))]
## tmp2$O.prev <- O.fwd[,Pm.prior[[j]][[prev]][[k]]$linked,,]
## tmp2$prior <- Pm.prior[[j]][[prev]][[k]]
## tmp2$linked <- length(Pm.prior[[j]][[prev]][[k]]$linked) > 0
## tmp2$P <- P[j,k]

tmp2 <- new.env()
tmp2$O.prev <- O.yri
tmp2$P <- Pm.prior[[j]]$freq['YRI']
tmp2$prior <- Pm.prior[[j]]$upstream$YRI
tmp2$linked <- length(linked.yri) > 0

##### Set up environment for P.gammas.a #####
tmp <- new.env()
tmp$fun <- 'P.gammas.a'
tmp$geno <- haps[,j,]
tmp$haps <- haps
tmp$prior <- Pm.prior[[j]]$model
tmp$P <- P[j,]
tmp$lambda <- lambda
tmp$d <- d[j]
## tmp$d <- d[j+1]
tmp$Ak <- Ak
tmp$gammas.prev <- gammas.prev

plot(ps[,,1], Ak[,,1]); abline(0:1)
plot(ps[,,2], Ak[,,2]); abline(0:1)

plot(pag[,,1]*pg[,,1], pag[,,2]*pg[,,2])

plot(ps[,,1], globalenv()$gammas[,j,,1])

plot(Ak[,,1], as.matrix(anc[,c('A1', 'A2')]), col = rep(c('red', 'blue'), each = dim(anc)[1])); abline(0:1)
plot(Ak[,1,1], Ak[,2,1])

##### Set up environment for updt.A #####

tmp2 <- new.env()
tmp2$A <- clusterEvalQ(tmp$cl, A0)[[1]]
tmp2$Ak <- clusterEvalQ(tmp$cl, Ak)[[1]]
tmp2$Ajs <- clusterEvalQ(tmp$cl, gammas[, chr != sex.chr,,])[[1]]
tmp2$omega <- clusterEvalQ(tmp$cl, omega)[[1]]
tmp2$sigma <- clusterEvalQ(tmp$cl, sigmaA)[[1]]
tmp2$epsilon <- 0.01
tmp2$dev <- clusterEvalQ(tmp$cl, dev)[[1]]

tmp <- new.env()
tmp$fun <- 'updt.A'
tmp$A <- A0
tmp$Ak <- Ak
tmp$G <- G[, chr != sex.chr,,]
tmp$omega <- omega
tmp$sigma <- sigmaA
tmp$epsilon <- 0.01
tmp$phased <- !is.null(haps)
tmp$dev <- dev

##### Set up environment for updt.lambda #####
tmp <- new.env()
tmp$fun <- 'updt.lambda'
tmp$G <- clusterEvalQ(cl, G[,chr != sex.chr,,])[[1]]
tmp$G0 <- clusterEvalQ(cl, G0[,chr != sex.chr,])[[1]]
tmp$A0 <- clusterEvalQ(cl, A0)[[1]]
tmp$Ak <- clusterEvalQ(cl, Ak)[[1]]
tmp$lambda <- clusterEvalQ(cl, lambda)[[1]]
tmp$sigma <- clusterEvalQ(cl, sigmaL)[[1]]
tmp$distances <- clusterEvalQ(cl, d[chr != sex.chr])[[1]]
tmp$alpha <- clusterEvalQ(cl, alpha)[[1]]
tmp$beta <- clusterEvalQ(cl, beta)[[1]]
tmp$chr <- clusterEvalQ(cl, chr[chr != sex.chr])[[1]]
tmp$P.Aeq1 <- clusterEvalQ(cl, P.Aeq1)[[1]]
tmp$dev <- clusterEvalQ(cl, dev)[[1]]
tmp$W = NULL

tmp <- new.env()
tmp$fun <- 'updt.lambda'
tmp$G <- G[,chr != sex.chr,,, drop = FALSE]
tmp$G0 <- G0[,chr != sex.chr,]
tmp$A0 <- A0
tmp$Ak <- Ak
tmp$lambda <- lambda
tmp$sigma <- sigmaL
tmp$distances <- d[chr != sex.chr]
tmp$alpha <- alpha[1]
tmp$beta <- alpha[2]
tmp$chr <- chr[chr != sex.chr]
tmp$P.Aeq1 <- P.Aeq1
tmp$phased <- !is.null(haps)
tmp$dev <- dev
tmp$W = NULL

# for X chromosome
tmp2 <- new.env()
tmp2$G <- clusterEvalQ(cl, G[,chr == sex.chr,,])[[1]]
tmp2$G0 <- clusterEvalQ(cl, G0[,chr == sex.chr,])[[1]]
tmp2$A0 <- clusterEvalQ(cl, AX0)[[1]]
tmp2$Ak <- clusterEvalQ(cl, AX)[[1]]
tmp2$lambda <- clusterEvalQ(cl, lambdaX)[[1]]
tmp2$sigma <- clusterEvalQ(cl, sigmaL)[[1]]
tmp2$distances <- clusterEvalQ(cl, d[chr == sex.chr])[[1]]
tmp2$alpha <- clusterEvalQ(cl, alpha)[[1]]
tmp2$beta <- clusterEvalQ(cl, beta)[[1]]
tmp2$chr <- clusterEvalQ(cl, chr[chr == sex.chr])[[1]]
tmp2$P.Aeq1 <- clusterEvalQ(cl, P.Aeq1)[[1]]
tmp2$dev <- clusterEvalQ(cl, dev)[[1]]
tmp2$W <- clusterEvalQ(cl, W)[[1]]


##### Set up environment for logP.delta.sigma #####
tmp3 <- new.env()
tmp3$deltas <- tmp2$deltas.new
tmp3$sigma <- tmp2$sigma
tmp3$lims <- tmp2$lims
tmp3$dev <- tmp2$dev
tmp3$epsilon <- 0.01


##### Set up environment for updt.sigma #####
tmp3 <- new.env()
tmp3$sigma <- tmp2$sigma
tmp3$deltas <- tmp2$deltas
tmp3$lims <- tmp2$lims
tmp3$dev <- tmp2$dev
tmp3$burn <- 50
tmp3$iter <- 100
tmp3$isA <- FALSE


##### Set up environment for sample.lambda #####
tmp <- new.env()
tmp$fun <- 'sample.lambda'
tmp$k <- k
tmp$G <- G[include,,,, drop = FALSE]
tmp$Ak <- Ak[include,,]
tmp$lambda <- lambda[include,]
tmp$Pm.prior <- Pm.prior
tmp$distances <- distances
tmp$chr <- chr
tmp$alpha <- alpha
tmp$beta <- beta
## tmp$P.Aeq1 <- P.Aeq1
tmp$phased <- phased
tmp$dev <- dev



##### short simulation of lambda and A #####

individuals <- data.frame(Afri = c(rep(1, 82), rep(0, 18)),
                          Euro = c(rep(0, 82), rep(1, 18)),
                          lambda = rep(0, 100),
                          count = rep(1, 100))

par(mar = c(5, 5, 2, 2))
plot(rep(0, 2), 0:1, xlab = expression(lambda), ylab = expression(A[Afri]), xlim = c(0, 10), ylim = 0:1,
     cex.lab = 2, cex.axis = 1.5, bty = 'l', pch = 20, col = rgb(0,0,0,.2))

while(!'stop' %in% system('ls', intern = TRUE))
{
    who <- sample(c(1:dim(individuals)[1])[individuals$count == 1], size = 2)
    individuals$count[who] <- 0

    new.euro <- mean(individuals$Euro[who])
    new.afri <- mean(individuals$Afri[who])

    new.lambda <- with(individuals, mean(ifelse(Euro[who] == 1 | Afri[who] == 1, 0, lambda[who] + 1)))

    individuals <- rbind(individuals, c(new.afri, new.euro, new.lambda, 1))
    points(new.lambda, new.afri, pch = 20, col = rgb(0,0,0,.2))
}


##### Figures of stuff from cluster #####

par(mfrow = c(2, 4))

invisible(lapply(clusterEvalQ(cl, A0), function(x) hist(x[,1])))
invisible(lapply(clusterEvalQ(cl, lambda), function(x) hist(x[,1])))


##### updt.tau.Pm ###
tmp <- new.env()
tmp$fun <- 'updt.tau.Pm'
tmp$tau <- tau
tmp$Pm <- Pm
tmp$Pm.counts <- Pm.counts
tmp$Pm.prior <- Pm.prior
tmp$P <- P
tmp$to.update <- to.update
tmp$burn <- 50
tmp$iter <- 100
tmp$every <- 10
tmp$debug <- FALSE
tmp$dev <- dev

##### sample.tau.Pm.prior #####
tmp <- new.env()
tmp$fun <- 'sample.tau.Pm.prior'
tmp$prior <- Pm.prior[[j]]$model[[k]]
tmp$tau <- tau[k]
tmp$tau.new <- tau.new
tmp$ev <- to.update[[k]]$priors[[jj]]$ev
tmp$dev <- dev

##### updt.sigma2 #####
tmp <- new.env()
tmp$fun <- 'updt.sigma2'
tmp$prior <- prior
tmp$tau <- tau
tmp$loglik.cur <- loglik.cur
tmp$burn <- 50
tmp$iter <- 100
tmp$every <- 10

##### initialize.G #####
tmp <- new.env()
tmp$fun <- 'initialize.G'
tmp$Pm.prior <- Pm.prior
tmp$lambda <- lambda
tmp$chr <- chr
tmp$Ak <- Ak
tmp$pos <- pos
