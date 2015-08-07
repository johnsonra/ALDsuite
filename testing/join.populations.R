# join.populations.R
# Join disperate files created in make.populations.R into one data file
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created May 6, 2014
# Last Modified May 9, 2014

library(rhdf5)

# get a list of all files to be merged
files <- system('ls pops | grep .h5', intern = TRUE)

groups <- gsub('.h5', '', gsub('Population', '', files), fixed = TRUE)

files <- paste('pops/', files, sep = '')

# repository for anc for everyone
anc <- NULL

# For each group:
for(i in 1:length(groups))
{
    # get anc
    tmp <- h5read(files[i], 'anc')
    colnames(tmp) <- c('id', 'group', 'A1', 'A2', 'A0', 'male')

    if(is.null(anc))
    {
        anc <- tmp
    }else{
        anc <- rbind(anc, tmp)
    }
}

anc <- as.data.frame(anc)

save(anc, file = 'pops/anc.RData')
