# gen.calc.R
# Calculate genetic position based on Matise et al map
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 14, 2006
# Last Modified October 30, 2013

gen.calc <- function(chrom, pos, n.extrap = 20, map = 'rutgers37', warn = FALSE)
{
    if(length(chrom) != length(pos))
        stop("chrom and pos must be of the same length.")

    if(n.extrap > 50 | n.extrap < 3)
        stop("n.extrap must be between 3 and 50.")

    # load appropriate rutgers map or throw an error if the map doesn't exist
    switch(map,
           rutgers37 = {
                            data(rutgers37)
                            rutgers <- subset(rutgers37,
                                              select = c(phys.pos, cM, chr))
                            names(rutgers)[2] <- 'gen.pos'
                       },
           rutgers36 = {
                            data(rutgers36)
                            rutgers <- subset(rutgers36,
                                              select = c(phys.pos, cM, chr))
                            names(rutgers)[2] <- 'gen.pos'
                       },
           stop(map, ' map not currently supported'))

    # make sure things are ordered in a good way...
    # this is neccesary due to the way we compile the results
    # use rev.order to undo this in the end
    chrom.order <- order(chrom, pos)
    rev.order <- order((1:length(chrom))[chrom.order])

    chrom <- chrom[chrom.order]
    pos <- pos[chrom.order]

    # get neighbors and interpolate genetic position (by chromosome)
    gen.pos <- NULL

    for(i in unique(chrom[chrom %in% 1:23]))
    {
      rutgers.sub <- rutgers[rutgers$chr == i,]

      # get positions for this chromosome
      gen.pos.sub <- unlist(sapply(pos[chrom == i], function(x)
                     {
                       # find distance to all neighbors
                       rutgers.sub$dist <- rutgers.sub$phys.pos - x

                       # check if the position is given in rutgers map
                       # I take the mean here, because 4 pairs of markers
                       # (micro satellites) end up mapped to the same
                       # physical postion, but different genetic positions
                       if(any(rutgers.sub$dist == 0))
                         return(with(rutgers.sub, mean(gen.pos[dist == 0])))

                       # check that there is a neighbor on both sides
                       if(with(rutgers.sub, min(dist) < 0 & max(dist) > 0))
                       {
                         # get neigbors to either side
                         rutgers.sub <- neighbors(rutgers.sub)

                         # calculate slope
                         m <- with(rutgers.sub, (gen.pos[2] - gen.pos[1]) /
                                                (phys.pos[2] - phys.pos[1])
                                  )

                         # return estimate by point slope formula
                         return(with(rutgers.sub, m*(x - phys.pos[1]) + gen.pos[1]))
                       }else{
                         if(warn)
                           warning(paste('Extrapolating genetic map position for',
                                         'chromosome', i, 'at position', x))

                         ### makes our values non-decreaseing and force of the intercept through 0 ###
                         if(all(rutgers.sub$dist > 0))
                         {
                             # force minimum point to the origin at the beginning of the chromosome
                             phys.mod <- min(rutgers.sub$phys.pos)
                             gen.mod <- min(rutgers.sub$gen.pos)
                         }else{

                             # force maximum point to the origin at the end of the chromosome
                             phys.mod <- max(rutgers.sub$phys.pos)
                             gen.mod <- max(rutgers.sub$gen.pos)
                         }

                         rutgers.sub$phys.pos <- rutgers.sub$phys.pos - phys.mod
                         rutgers.sub$gen.pos <- rutgers.sub$gen.pos - gen.mod

                         # get nearest n neighbors
                         rutgers.sub <- neighbors2(rutgers.sub, n.extrap)

                         # estimate line
                         rutgers.pred <- lm(gen.pos ~ 0 + phys.pos, data=rutgers.sub)

                         return(predict(rutgers.pred, newdata=data.frame(phys.pos=x - phys.mod)) + gen.mod)
                       }
                     }))

      # compile results
      gen.pos <- c(gen.pos, gen.pos.sub)
    }

    # return results
    return(data.frame(chr = chrom[rev.order], phys.pos = pos[rev.order], gen.pos = gen.pos[rev.order]))
}

# find nearest 2 neigbors on either side of a point, distance to the point is contained in dat
neighbors <- function(dat)
{
    # define groups on either side of 0
    dat$gt0 <- dat$dist > 0
    dat$lt0 <- dat$dist < 0

    # order the groups so we can identify the one closest to 0.
    # note that we have a few microsatelites mapped to the same
    # physical position, but at different genetic positions that
    # we need to watch out for.
    gt0.order <- with(dat, order(dist[gt0], gen.pos[gt0]))
    lt0.order <- with(dat, order(-dist[lt0], -gen.pos[lt0]))

    # return physical and genetic positions for closest to 0
    dat$keep <- FALSE
    dat$keep[dat$gt0][gt0.order[1]] <- TRUE
    dat$keep[dat$lt0][lt0.order[1]] <- TRUE

    return(subset(dat, keep, select=c('phys.pos', 'gen.pos')))
}

# find nearest 4 neighbors to a point, distance to the point is contained in dat
neighbors2 <- function(dat, n.extrap)
{
    # order data so we can identify the ones closest to 0.
    # note that in this case 0 is assumed to be a bound of dist.
    # note also that we have a few microsatelites mapped to the
    # same physical position, but at different genetic positions
    # that we need to watch out for.
    dat.order <- with(dat, order(abs(dist), sign(dist) * gen.pos))

    # return physical and genetic positions for closest n
    dat$keep <- FALSE
    dat$keep[dat.order[1:n.extrap]] <- TRUE

    return(subset(dat, keep, select=c('phys.pos', 'gen.pos')))
}
