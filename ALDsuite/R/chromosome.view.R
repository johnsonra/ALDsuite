# chromosome.view.R
# Plot the chromosome view of an amap object
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC-Frederick, Inc
# Created August, 2005
# Last Modified April 4, 2008


chromosome.view <- function(x, chr.len, note, add.legend, xaxt, type, bty,
                            main, ...)
{
  # set default graphic parameters if not already specified
  if(is.null(xaxt))
    xaxt <- 'n'
  if(is.null(bty))
    bty <- 'n'

##################
# Ethnicity Data #
##################

# make sure there is data to plot
  if(!is.null(x$stats$marker.stats))
  {
    # we only want to map the actual markers
    marker.stats <- subset(x$stats$marker.stats,
                           substr(as.character(x$stats$marker.stats$snp.id), 1, 4) != 'fake')

    marker.stats$phys.pos <- marker.stats$phys.pos

    # Plot ethnicity data
    theta.mean <- x$stats$estimates$theta.mean
    thetax.mean <- x$stats$estimates$thetax.mean

    marker.stats$line1 <- (marker.stats$g.case -
                           ifelse(marker.stats$chr == 23, thetax.mean,
                                  theta.mean)) * 100
    marker.stats$line2 <- (marker.stats$g.control -
                           ifelse(marker.stats$chr == 23, thetax.mean,
                                  theta.mean)) * 100

    # should we be plotting african ancestry here? that is what the label says.
    if(theta.mean < .5)
    {
      marker.stats$line1 <- -marker.stats$line1
      marker.stats$line2 <- -marker.stats$line2
    }

    teal <- rgb(19, 132, 130, maxColorValue = 255)
    .plot.amap2(marker.stats, to.plot = c('case', 'control'),
                color = c(rgb(0,0,1), teal), note, eth=TRUE,
                add.legend = add.legend, xaxt = xaxt, type = type, bty = bty,
                main, x$genome.score, chr.len, x$stats$chr.stats, ...)
  }

########################
# lgd / ccs Statistics #
########################

##   if(!is.null(stats$marker.stats))
##   {
##     marker.stats <- stats$marker.stats
##     marker.stats$phys.pos <- NULL # we only want to map the actual markers

##     marker.stats <- merge(marker.stats, subset(markers,
##                                                select=c('snp.id', 'phys.pos')))
##     ylim <- c(min(c(marker.stats$ccs, marker.stats$lgs), na.rm=TRUE),
##               max(c(marker.stats$ccs, marker.stats$lgs), na.rm=TRUE))
##     marker.stats$line1 <- marker.stats$lgs
##     marker.stats$line2 <- marker.stats$ccs

##     .plot.amap2(marker.stats, main, ylim, ylab='LOD Score',
##                to.plot=c('lgs', 'ccs'), color=c('green', 'orange'),
##                genome.score)
##   }

}

.plot.amap2 <- function(dat, to.plot, color, note, eth, add.legend, xaxt, type,
                        bty, main, genome.score, chr.len, chr.stats=NULL, ...)
{
  beg <- .07 # start the plots here (in [0,1])

  ylim <- range(c(dat$line1, dat$line2))

# plot titles etc...
  par(fig=c(0,1,0,1), mar=c(2,2,4,0))
  plot(NA, NA, xlab='', ylab='', bty=bty, xaxt=xaxt, yaxt='n', xlim=c(0,1),
       ylim=c(0,1))

  if(!is.null(main))
    mtext(main, cex=2.5, line=2)

  if(eth){
    mtext(c('Change in African', 'ancestry vs. genome average (%)'),
          side = 2, line = c(3,2), las = 0)
  }else{
    mtext('LOD Score', side=2, line=.5, cex=1.5)
  }
  mtext('Position in Mb', side=1, line=1, cex=1.5)

  if(add.legend)
  {
    legend(.8, .2, legend=to.plot, col=color, lty=1, bty='n', cex=2)
    text(.89, .37, paste('risk: ', genome.score$risk1[1], sep=''), adj=0,
         cex=1.5)
    text(.89, .3, paste('lgs: ', genome.score$score[1], sep=''), adj=0, cex=1.5)
  }

  par(mar = c(2.6,0,0,0))

##   # add note
##   if(!is.null(note))
##     text(xlim[1], ylim[1] + ycorrect, note, pos=4)


# plot each chromosome
  for(i in 1:23)
  {
    if(i %in% c(1, 5, 10, 16))
    {
      # set plot boundaries for each line
      if(i == 1){
        y <- c(.715, .94)
      }else{
      if(i == 5){
        y <- c(.49, .715)
      }else{
      if(i == 10){
        y <- c(.265, .49)
      }else{
        y <- c(.04, .265)
      }}}

      # reset end to beg
      end <- beg

      # for use in setting up yaxis
      yaxt <- 's'
    }else{
      yaxt <- 'n'
    }

    tmp <- subset(dat, dat$chr == i)
    x <- c(end, chr.len[i] / 1000 + end - .01)
    end <- chr.len[i] / 1000 + .01 + end # save our place for next chromosome
    par(fig = c(x, y), new=TRUE)

    # plot chromosome
    plot(tmp$phys.pos[order(tmp$phys.pos)], tmp$line1[order(tmp$phys.pos)], type='l',
         col=color[1], bty='n', ylim=ylim, yaxt=yaxt, xaxt='n',
         xlim=c(0, chr.len[i]))

    # label x axis
    axis(1, tick=FALSE)

    # show snp density
    axis(1, at=tmp$phys.pos[order(tmp$phys.pos)], labels=FALSE)

    # if we have data for two lines, plot the second
    if(length(to.plot) > 1)
      lines(tmp$phys.pos[order(tmp$phys.pos)], tmp$line2[order(tmp$phys.pos)], col=color[2])

    # label chromosome
    mtext(paste('Chr', i, sep=' '), side=1, line=1.7, cex=.8)

    # Add chromosome statistics
    if(!is.null(chr.stats))
    {
      text(0, .8, label=chr.stats$lgs.max[i], adj=0, cex=.75)
    }else{
      if(i %in% c(4,9,15,23)) # Add significance line for marker stats
      {
        par(fig = c(beg, end + .01, y), new=TRUE)
        plot(c(-1,1), c(2,2), xlim=c(0,1), ylim=ylim, xaxt='n', yaxt='n',
             type='l', col='grey75', bty='n')
      }
    }
  }
}
