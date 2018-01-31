# show.chain.R
# Show the state of a variable in the MCMC chain over the course of an analysis
# Randall Johnson
# CCR Collaborative Bioniformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc

# this function takes the following arguments:
#     ncpath - path to a NetCDF file created by admixture()
#     param  - character name of the parameter to be shown
#     ...    - any other graphical parameters the user wishes to have passed to plot
show.chain <- function(ncpath, param, ...)
{

    ######### Checks #########
    # check that we have the correct library to access the nc file
    if(!require(ncdf4, quietly = TRUE))
        stop("ncdf4 package required for this function")

    # check that the path provided is a good one
    if(!file.exists(ncpath))
        stop("There seems to be a problem with the provided path")


    ######### Setup #########
    # open file and get the parameter in question
    debug <- nc_open(ncpath)

    chain <- ncvar_get(debug, param)

    # the last dimension is for itterations of the MCMC chain -- find out how many dimensions we are working with
    dims <- dim(chain)

    # figure out how many sparklines to plot on each page
    nfigs <- prod(dims[-length(dims)])

    if(nfigs <= 5)
    {
        ncol <- 1
        nrow <- nfigs
    }

    if(nfigs > 5 & nfigs <= 20)
    {
        ncol <- 1
        nrow <- 5
    }

    if(nfigs > 20 & nfigs <= 100)
    {
        ncol <- 3
        nrow <- 5
    }

    if(nfigs > 100)
    {
        ncol <- 10
        nrow <- 5
    }

    # add a few defaults to plot that aren't a part of par()
    defaults <- "..."

    if(!'type' %in% names(list(...)))
        defaults <- c(defaults, "type = 'b'")

    if(!'ylab' %in% names(list(...)))
        defaults <- c(defaults, "ylab = param")


    ######### Plotting #########
    # plot each line in the itteration dimension
    pdf(paste(param, '.pdf', sep = ''), width = 11, height = 8.5)
    par(cex = .5, bty = 'n', las = 1, mfrow = c(nrow, ncol))

    eval(parse(text = paste("apply(chain, 1:(length(dims) - 1), plot,",
                            paste(defaults, collapse = ', '), ")")))


    dev.off()
}
