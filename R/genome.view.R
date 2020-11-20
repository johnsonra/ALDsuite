# genome.view.R
# Show the genome-wide view of an amap object
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created April 2, 2008
# Last Modified August 14, 2012


genome.view <- function(x, cM, add.legend, change.par, chr.len, ...)
{
  args <- list(...)

  # check/set defaults
  if(is.null(args$xaxt))
    args$xaxt <- 'n'
  
  if(is.null(args$bty))
    args$bty <- 'n'

  if(is.null(args$type))
    args$type <- 'l'

  if(is.null(args$col))
    color <- c('blue', rgb(19, 132, 130, maxColorValue = 255),#teal
               'red', 'black')
  else
    color <- rep(args$col, 4)[1:4]

  if(is.null(args$xlab))
    args$xlab <- ifelse(cM, 'Genetic Position in cM', 'Position in Mb')

  if(is.null(args$ylab))
    ylab <- c('LOD Score', 'Change in African',
              'ancestry vs. genome average (%)')
  else
    ylab <- c(rep(args$ylab, 2)[2:1], '')
      
  with(x,
  {
    if(!is.null(stats$marker.stats))
    {
      if(change.par)
      {
        par.args <- par(no.readonly = TRUE) # save this to restore later
        par(mfcol=c(2,1), mar=c(5.1, 4.1, 0, 1.1), las=2)
      }

      # we only want to map the actual markers
      marker.stats <- subset(stats$marker.stats,
                             substr(as.character(snp.id), 1, 4) != 'fake')

      marker.stats$phys.all <- marker.stats$phys.pos / 10^6 +
        sapply(marker.stats$chr, function(x){sum(chr.len[0:(x-1)])})

      # Plot ethnicity data
      theta.mean <- stats$estimates$theta.mean
      thetax.mean <- stats$estimates$thetax.mean

      marker.stats$line1 <- (marker.stats$g.case -
                             ifelse(marker.stats$chr == 23, thetax.mean,
                                    theta.mean)) * 100
      if(all(marker.stats$g.control < 0.001)){
        marker.stats$line2 <- NA
      }else{
        marker.stats$line2 <- (marker.stats$g.control -
                               ifelse(marker.stats$chr == 23, thetax.mean,
                                      theta.mean)) * 100
      }

      # should we be plotting african ancestry here?
      # That is what the y-label says
      if(theta.mean < .5)
      {
        marker.stats$line1 <- -marker.stats$line1
        marker.stats$line2 <- -marker.stats$line2
      }

      # run plotting subfunction for case/control frequencies
      args$dat <- marker.stats
      if(all(is.na(marker.stats$g.control))){# double check that case/control data exists...
        args$to.plot <- 'case'
        args$colors <- color[1]
      }else{
        args$to.plot <- c('case', 'control')
        args$colors <- color[1:2]
      }
      args$eth <- TRUE
      args$add.legend <- add.legend
      args$chr.len <- chr.len
      args$ylab <- ylab[2:3]
      do.call('.plot.amap', args)

      # Plot statistics
      # get the associated p from the z score
      tmp.p <- 2*pnorm(abs(marker.stats$ccs), lower.tail = FALSE)
      # translate it into a LOD score
      marker.stats$line2 <- -log10(tmp.p / (1-tmp.p)) # this will be NA if missing ccs (eg it is 0)
      marker.stats$line1 <- marker.stats$lgs

      if(all(is.na(marker.stats$line2))){# Again, does the ccs exist???
        args$to.plot <- 'case'
        args$colors <- color[3]
      }else{
        args$to.plot <- c('case', 'control')
        args$colors <- color[3:4]
      }

      # run plotting subfunction for statistics
      args$dat <- marker.stats
      args$to.plot <- c('lgs', 'ccs')
      args$eth <- FALSE
      args$ylab <- c(ylab[1], '')
      do.call('.plot.amap', args)
    }
  
    if(change.par) # reset par
      lapply(names(par.args), function(x)
             eval(parse(text = paste('par(', x, '=',
                          paste(deparse(par.args[[x]]), collapse = ''), ')'))))

  })# end with(x, ...
}

.plot.amap <- function(dat, to.plot, eth, add.legend, chr.len, colors, ...)
{
  # start point of each chromosome (and end of 23)
  chr.sum <- sapply(1:24, function(x){sum(chr.len[0:(x-1)])})

  args <- list(...)
  ylab <- args$ylab # save for later...

  if(!is.null(args$sub))
  {
    args$note <- args$sub
    args$sub <- NULL
  }else{
    args$note <- NULL
  }
    
  if(is.null(args$ylim))
    args$ylim <- with(dat, c(floor(min(c(line1[line1 > -Inf], line2[line2 > -Inf]),
                                       na.rm=TRUE)),
                             ceiling(max(c(line1[line1 < Inf], line2[line2 > Inf]),
                                         na.rm=TRUE))))
  
  yrange <- args$ylim[2] - args$ylim[1]
  
  args$lab <- c(40, max(round(yrange) / 2^floor(yrange/10), 3), 7)
    
  # xlim not settable by end user
  xlim <- range(dat$phys.all)
  # trim xlim to get rid of extra whitespace on either end of plot
  xlen <- xlim[2] - xlim[1]

  trim.factor <- 1 - 1/1.038 # 4% is added, this is about what I want to trim
  xlim <- xlim + trim.factor*xlen*c(1,-1)

  # plot data and label axes
  args$x <- dat$phys.all
  args$y <- dat$line2
  args$col <- colors[2]
  args$ylab <- ylab[1]
  do.call('plot', args)

##     if(!is.null(args$h)) #this is messed up...default works, but not settable
##       abline(h = args$h, col='grey75')
##     else
##       abline(h = 0, col = 'grey75')

  args$y <- dat$line1
  args$col <- colors[1]
  do.call('lines', args)

  # add second line of ylabel if there is one
  mtext(ylab[2], side=2, line=2, las=0)

  # set up x axis
##     xaxis <- data.frame(points = axis(1, labels=FALSE, tick=FALSE, lty=0))

##     xaxis$chr <- as.numeric(cut(xaxis$points, breaks=chr.sum))

##     xaxis$chr.pos <- paste(xaxis$chr, round(xaxis$points - chr.sum[xaxis$chr]),
##                            sep='-')

##     # print x axis -- need to add some here!
##     axis(side=1, at=xaxis$points, labels=xaxis$chr.pos, lty=0, tick=FALSE,
##          line=-.4)
  axis(side=1, at=dat$phys.all, tick=TRUE, labels=FALSE)
    
  # print chromosomes
  ycorrect <- yrange/40
  ymin <- par()$usr[3] + ycorrect
  xrange <- par()$usr[2] - par()$usr[1]
  
  xcorrect <- c(1, -1)*xrange/250 # to make up for rounded end
    
  for(i in 1:23)
    lines(chr.sum[i:(i+1)] + xcorrect, rep(ymin, 2), lwd=7, col='grey50',
          lend = 'round')

  # add note
  if(!is.null(args$note))
    text(xlim[1], args$ylim[1] + ycorrect, args$note, pos=4)

  # add legend
  if(add.legend)
    legend('topleft', to.plot, col=colors, lty=1, bty='n')
}
