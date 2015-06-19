# data.R
# Creation of example datasets for inclusion in ALDsuite
# Randall Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc


library(devtools)

source_url('https://raw.githubusercontent.com/johnsonra/ALDdata/dev/R/helpers.R')

# read in and save chromosome 22 for ASW cohort
asw <- get.phased(22, 'ASW', dataDir = '../ALDdata/Data/')

save(asw, file = 'ALDsuite/data/asw.RData')
