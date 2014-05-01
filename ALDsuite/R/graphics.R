# graphics.R
# Graphics functions for ALDsuite objects
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created April 30, 2014
# Last Modified April 30, 2014


# x = gammas object from admixture object

ALDplot.gammas <- function(x, indiv = NULL, phased = TRUE, col = NULL)
{
    if(!phased)
    {
        warning('Unphased plotting not yet implemented')
        return()
    }

    # if no individuals are specified, plot all of them
    if(is.null(indiv))
        indiv <- dim(x)[1]

    # default colors
    if(is.null(col))
    {
        if(phased)
        {
            col <- 1:dim(x)[4]
        }else{
            stop('still need default colors')
        }
    }

    for(i in indiv)
    {
        moms <- NULL
        dads <- NULL

        for(l in 1:dim(x)[4])
        {
            moms <- cbind(moms, x[i,,1,l])
            dads <- cbind(dads, x[i,,2,l])
        }

        barplot(rbind(t(dads), t(moms)), space = 0, axisnames = FALSE, axes = FALSE,
                col = rep(col, 2), border = NA)
    }
}
