# tools.R
# General tools for use in ALDsuite
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created September 3, 2013
# Last Modified September 3, 2013

# tested...it works! :)
categorize <- function(mat, ref, long = TRUE, collapse = TRUE)
{
    if(class(mat) == 'data.frame')
        mat <- as.matrix(mat)

    # this will make the code nicer below
    if(!long)
        mat <- t(mat)

    # create the new matrix (numeric)
    newmat <- matrix(0, nrow = dim(mat)[1], ncol = dim(mat)[2], dimnames = list(rownames(mat), colnames(mat)))

    # populate newmat
    for(i in 1:dim(mat)[1])
    {
        for(j in 1:dim(mat)[2])
        {
            if(mat[i,j] != ref[i])
                newmat[i,j] <- 1
        }
    }

    # if we are collapsing, do so now
    if(collapse)
    {
        collapsed <- matrix(nrow = dim(mat)[1], ncol = dim(mat)[2] / 2,
                            dimnames = list(rownames(mat), colnames(mat)[1:(dim(mat)[2] / 2) * 2 - 1]))

        for(i in 1:dim(collapsed)[1])
        {
            for(j in 1:dim(collapsed)[2])
            {
                collapsed[i,j] <- newmat[i,j*2 - 1] + newmat[i, j*2]
            }
        }

        if(!long)
            return(t(collapsed))

        return(collapsed)
    }

    if(!long)
        return(t(newmat))

    return(newmat)
}

# C implementation of the above

catagorize.C <- function(mat, ref, long = TRUE, collapse = TRUE)
{
    categorize(mat, ref, long, collapse)
}
