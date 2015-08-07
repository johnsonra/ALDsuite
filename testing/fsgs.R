# fsgs.R
# Run of FSGS data
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created December 13, 2013
# Last Modified December 13, 2013

library(ALDdata)
library(ALDsuite)

# read in genotypes
genos <- unique(read.csv('FSGSgenotypes.txt', stringsAsFactors = FALSE))
clin <- unique(subset(genos, select = c(hgal, race, sex, fsgs, hiv)))
genos <- subset(genos, select = c(hgal, db.snp.no, genotype))

genos <- reshape(genos, v.names = 'genotype', timevar = 'db.snp.no', idvar = 'hgal', direction = 'wide')
names(genos) <- gsub('genotype.', '', names(genos), fixed = TRUE)

genos <- as.matrix(genos)
rownames(genos) <- genos[,1]
genos <- genos[,-1]

# split into quasi haplotypes
dups <- matrix('', nrow = dim(genos)[1] * 2, ncol = dim(genos)[2], dimnames = list(NULL, colnames(genos)))

dads <- 1:dim(genos)[1] * 2
moms <- dads - 1

for(j in 1:dim(genos)[2])
{
    dups[moms,j] <- substr(genos[,j], 1, 1)
    dups[dads,j] <- substr(genos[,j], 2, 2)
}

# double check conversion
data(hapmap)
hapmap <- subset(hapmap, rs %in% colnames(dups))

dups <- dups[,hapmap$

# convert into 0,1,2s
