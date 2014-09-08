# ibd.R
# Functions to infer relatedness
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Last Modified September 8, 2014

# formulas in this file follow those described in PLINK paper (Purcell 2007)

############################################################
# Observed number of shared alleles at the given IBS level #
############################################################

O.shared <- function(genos, ibs)
{
    # number of alleles shared
    shared <- 2 - abs(genos[1,] - genos[2,])

    return(sum(shared == ibs))
}


###################################################################
# IBD State across the entire genome for each pair of individuals #
###################################################################

# probably don't need to run this at every marker...use a subset (one per window)
# -- will need to be sure they are sufficiently frequent, though

# move this to C at some point...could be a lot faster and a little more memory efficient
ibd.pairs <- function(geno)
{
    # tally number of non-missing genotypes
    n <- apply(!is.na(geno), 2, sum)

    # tally number of alleles
    X <- apply(geno, 2, sum, na.rm = TRUE)
    Y <- 2*n - X

    # avoid negative numbers (these will be nearly non-informative anyway, so just skip them)
    todrop <- which(X < 3 | Y < 3)

    n <- n[-todrop]
    X <- X[-todrop]
    Y <- Y[-todrop]

    # allele frequencies
    lp <- log(X / (2*n))
    lq <- log(Y / (2*n))

    # save some time on these where we need to calculate more than once
    logn <- log(n)
    logn_1 <- log(n - 1)
    logn_2 <- log(n - 2)
    logn_3 <- log(n - 3)

    lX <- log(X)
    lX_1 <- log(X - 1)
    lX_2 <- log(X - 2)

    lY <- log(Y)
    lY_1 <- log(Y - 1)
    lY_2 <- log(Y - 2)

# see Table 1 of Purcell 2007

    ######### E(IBS = 0 | IBD = 0) #########
    ibs0ibd0 <- sum(exp(2*lp + 2*lq + lX_1 - lX + lY_1 - lY + 3*logn - logn_1 - logn_2 - logn_3))

    ######### E(IBS = 1 | IBD = 0) #########
    ibs1ibd0 <- sum(exp(log(4) + 3*lp + lq + lX_1 + lX_2 - 2*lX + 3*logn - logn_1 - logn_2 - logn_3) +
                    exp(log(4) + lp + 3*lq + lY_1 + lY_2 - 2*lY + 3*logn - logn_1 - logn_2 - logn_3))

    ######### E(IBS = 2 | IBD = 0) #########
    ibs2ibd0 <- sum(exp(4*lp + lX_1 + lX_2 + lX_3 - 3*lX + 3*logn - logn_1 - logn_2 - logn_3) +
                    exp(4*lq + lY_1 + lY_2 + lY_3 - 3*lY + 3*logn - logn_1 - logn_2 - logn_3) +
                    exp(log(4) + 2*lp + 2*lq + lX_1 - lX + lY_1 - lY + 3*logn - logn_1 - logn_2 - logn_3))

    ######### E(IBS = 1 | IBD = 1) #########
    ibs1ibd1 <- sum(exp(log(2) + 2*lp + lq + lX_1 - lX + 2*logn - logn_1 - logn_2) +
                    exp(log(2) + lp + 2*lq + lY_1 - lY + 2*logn - logn_1 - logn_2))

    ######### E(IBS = 2 | IBD = 1) #########
    ibs2ibd1 <- sum(exp(3*lp + lX_1 + lX_2 - 2*lX + 2*logn - logn_1 - logn_2) +
                    exp(3*lq + lY_1 + lY_2 - 2*lY + 2*logn - logn_1 - logn_2) +
                    exp(2*lp + lq + lX_1 - lX + 2*logn - logn_1 - logn_2) +
                    exp(lp + 2*lq + lY_1 - lY + 2*logn - logn_1 - logn_2))

    ######### E(IBS = 2 | IBD = 2) #########
    ibs2ibd2 <- n
}
