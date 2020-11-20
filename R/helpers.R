# helpers.R
# A few helper functions for ALDsuite
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC Frederick, Inc
# Created September 10, 2013
# Last Modified September 26, 2013

# this will return a vector of one value per column (i.e. length(rows) should equal dim(mat)[2])
get.rows <- function(mat, rows)
{
    .Call("get_rows", as.numeric(mat), as.integer(rows - 1))
}

# tihs will return an approximate inverse for a matrix if it is not quite solvable
# code adapted from Hans W. Borchers in an Rhelp post
solve.approx <- function(X)
{
    Xinv <- try(solve(X), silent = TRUE)

    if(class(Xinv) == 'try-error')
    {
        s <- svd(X)
        Xinv <- s$v %*% diag(1/s$d) %*% t(s$u)
    }

    return(Xinv)
}


# this will reload the most recently loaded library
reload <- function()
{
    space <- search()[2]
    pack <- gsub('package:', '', space, fixed = TRUE)
    detach(unload = TRUE)
    eval(parse(text = paste('library(', pack, ')')))
}
