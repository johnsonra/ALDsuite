# ibd.R
# Functions to infer relatedness
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Last Modified September 6, 2014

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


