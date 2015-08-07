# prior.R
# Testing of setup.prior()
# Randall Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc

library(ALDdata)
library(ALDsuite)

# updates
updt <- parse(text = "{source('../ALDsuite/R/pcr.prior.R')
                       source('../ALDsuite/R/helpers.R')}")

# comparison/published code
comp <- parse(text = "")


##### One chromosome #####

data(asw)
markers <- asw$rsID[1:100]

eval(updt)
new <- setup.prior(markers, c('YRI', 'CEU'), unphased = FALSE)

eval(comp)
old <- setup.prior(markers, c('YRI', 'CEU'), unphased = FALSE)


##### Two chromosomes #####

data(yri20)
markers <- colnames(phased)[1:100]
data(yri21)
markers <- c(markers, colnames(phased)[1:100])

eval(updt)
new <- setup.prior(markers, c('YRI', 'CEU'), unphased = FALSE)
