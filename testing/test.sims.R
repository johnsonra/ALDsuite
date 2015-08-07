# test.sims.R
# Do all testing for tables
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created May 14, 2014
# Last Modified May 14, 2014


library(ALDsuite)


#########
# Setup #
#########

# get all files for this sample size
files <- system('ls pops | grep anc', intern = TRUE)
files <- files[files != 'anc.RData']

# strip out sample sizes (n)
tmp <- strsplit(files, split = '.', fixed = TRUE)
N <- sapply(tmp, length)

N <- ifelse(N == 4, sapply(tmp, `[`, 3),
                    sapply(tmp, `[`, 4))

# strip out seeds
tmp <- gsub('anc', '', files)

seeds <- sapply(strsplit(tmp, '.', fixed = TRUE), `[`, 1)

# stript out ORs
or <- tmp
for(i in 1:length(seeds))
    or[i] <- gsub(paste('.', N[i], '.RData', sep = ''), '', gsub(paste(seeds[i], '.', sep = ''), '', tmp[i]))

# add directory onto files
files <- paste('pops/', files, sep = '')

# power array
power <- array(0, dim = c(length(unique(or)), 3, 4),
               dimnames = list(unique(or), c(100, 200, 500), c(0.5, 0.75, 1.5, 2)))
power2 <- array(0, dim = c(length(unique(or)), 3),
                dimnames = list(unique(or), c(100, 200, 500)))


###########
# Testing #
###########

for(i in 1:length(files))
{
    load(files[i])

    if(!exists('out.phased'))
    {
        warning(paste('Skipping', i))
        next
    }
    ## for(phi in c(0.5, 0.75, 1.5, 2))
    ## {
    ##     lods <- lgs.case.only(out.phased, phi1 = phi, cases = anc$case, phased = TRUE)
    ##     power[or[i], '100', as.character(phi)] <- power[or[i], '100', as.character(phi)] + (max(lods) > 3) / 30
    ## }

    rownames(anc) <- anc$pid

    models <- mald(out.phased, anc)
    ps <- sapply(lapply(models, `[[`, 1), function(x)
                                              {
                                                  x <- try(x['g', 'Pr(>|z|)'])
                                                  if(class(x) == 'try-error')
                                                      return(1)
                                                  return(x)
                                              })

    power2[or[i], N[i]] <- power2[or[i], N[i]] + (min(ps, na.rm = TRUE) < 0.05 / length(models)) / 30

    rm(out.phased)
}
