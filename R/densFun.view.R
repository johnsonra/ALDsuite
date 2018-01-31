# plot.densFun.R
# plot the desnity function and 95% credible interval
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created May 22, 2007
# Last Modified August 8, 2008

densFun.view <- function(x, chromosome, ci, ci.col, ...)
{
  # use marker.stats for specified chromosome
  with(subset(x$stats$marker.stats, x$stats$marker.stats$chr == chromosome),
  {
    # set up plot
    args <- list(...)

    if(!is.null(args$xlim))
      args$xlim <- args$xlim * 1e6
    else
      args$xlim <- range(phys.pos)

    if(is.null(args$ylim))
      args$ylim <- range(10^lgs)

    if(is.null(args$main) & !is.na(ci))
      args$main <- paste(ci*100, "% Credible Interval", sep = '')
##     else if(is.null(args$main) & is.na(ci))
##       args$main <- paste("Relative Probability Density Function")

    if(is.null(args$xlab))
      args$xlab <- 'Physical Position (Mb)'

    if(is.null(args$ylab))
      args$ylab <- 'Locus-Genome Odds'

    args$x <- NA
    args$y <- NA
    args$xaxt <- 'n'

    do.call('plot', args)

    # calculate the probability density function
    dense.fun <- cum.prob.fun(gen.pos, 10^lgs)

    # 95% credible interval
    if(!is.na(ci))
    {
      if(!is.numeric(ci) | ci >= 1 | ci <= 0)
        stop("unless missing, ci must be numeric and in (0, 1)")

      from <- (1-ci)/2
      to <- 1 - from

      # get the intervals containing the end points
      ci.lower <- with(dense.fun, c(max(x[p <= from]), min(x[p >= from])))
      ci.upper <- with(dense.fun, c(max(x[p <= to]), min(x[p >= to])))

      # get the lower bound
      y <- unique(dense.fun$p[dense.fun$x %in% ci.lower])
      if(length(y) > 1) # make sure that we didn't hit a SNP dead on
      {
        m <- (y[1] - y[2]) / (ci.lower[1] - ci.lower[2])
        ci.lower <- ((from - y[1]) / m + ci.lower[1])*100 # if in morgans...
      }else{
        ci.lower <- y
      }

      # get the upper bound
      y <- unique(dense.fun$p[dense.fun$x %in% ci.upper])
      if(length(y) > 1) # make sure that we didn't hit a SNP dead on
      {
        m <- (y[1] - y[2]) / (ci.upper[1] - ci.upper[2])
        ci.upper <- ((to - y[1]) / m + ci.upper[1])*100 # if in morgans...
      }else{
        ci.upper <- y
      }

      # get both cM and Mb positions of 95% credible interval
      cred.int <- phys.calc(rep(chromosome, 2), c(ci.lower, ci.upper))
      rownames(cred.int) <- c('Lower', 'Upper')

      dens.lower <- dense.atX(phys.pos, 10^lgs, cred.int$phys.pos[1])
      dens.upper <- dense.atX(phys.pos, 10^lgs, cred.int$phys.pos[2])

      polygon(c(rep(cred.int$phys.pos[1], 2),
                phys.pos[phys.pos > cred.int$phys.pos[1] &
                         phys.pos < cred.int$phys.pos[2]],
                rep(cred.int$phys.pos[2], 2)),
              c(0, dens.lower,
                10^lgs[phys.pos > cred.int$phys.pos[1] &
                       phys.pos < cred.int$phys.pos[2]],
                dens.upper, 0), col = ci.col, border = NA)
    }else{
      cred.int <- NULL
    }

    # peak
    lines(phys.pos, 10^lgs)

    # real markers in data
    not.fake <- substr(snp.id, 1, 4) != 'fake'
    points(phys.pos[not.fake], 10^lgs[not.fake])

    args <- list(...)

    # put in xaxis?
    if(is.null(args$xaxt))
      args$xaxt <- 's'
    if(args$xaxt != 'n')
    {
      args$side <- 1
      args$labels <- FALSE
      args$at <- do.call('axis', args)

      args$labels <- args$at / 1e6
      do.call('axis', args)
    }

    # legend
    if(!is.na(ci))
      legend('topleft', legend = paste(c('Upper bound:', 'Lower bound:'),
                          round(cred.int$phys.pos[2:1] / 10e5, 1)),
             col = 'transparent', bty = 'n')

    invisible(cred.int)
  })
}

# calculate the density at a given position
dense.atX <- function(x, y, at)
{
  x1 <- max(x[x < at])
  x2 <- min(x[x > at])
  y1 <- y[x == x1]
  y2 <- y[x == x2]

  m <- (y2 - y1) / (x2 - x1)

  b <- y1 - m*x1

  return(m*at + b)
}

# given x and y, return the cumulative probability function
cum.prob.fun <- function(x, y)
{
  # checks
  if(length(x) != length(y))
    warning('length(x) != length(y)')

  # set up data frame
  dat <- na.omit(data.frame(x = x[order(x)], y = y[order(x)]))

  # this is where our integral is summed as we move along
  dens.fun <- 0

  # get the number of elements
  n <- length(dat$x)

  # if we only have 1 element, return 1 -- nonsensical...
  if(n < 2)
    return(1)

  # now integrate the polygonal representation
  for(i in 2:n)
  {
    y.range <- range(dat$y[(i-1):i])
    wdth <- dat$x[i] - dat$x[i-1]

    dens.fun[i] <- dens.fun[i-1] + wdth*y.range[1] +
                   wdth*(y.range[2] - y.range[1])/2
  }

  return(data.frame(x = dat$x,
                    p = dens.fun / dens.fun[n],
                    d = dens.fun))
}
