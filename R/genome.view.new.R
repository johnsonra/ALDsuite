# genome.view.R
# Updated genome.view function to plot genome-wide results from low density scans
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created August 29, 2012
# Last Modified August 29, 2012

## will take matrix for y (assumes chromosome values are the same for each column of y)
## when y is a matrix, y[,i] ~ x is plotted for all i in 1:dim(y)[2]
## plotting parameters will be used in iterative fashion (i.e. if length == dim(y)[2] each will be used once)

# x = chromosome positions to plot
# y = values to be plotted at each point along the chromosome (can be a matrix, in which case dim(y)[2] == length(x))
# chr = chromosome number for each x (assumed chromosomes are the same for all columns of x)
# conv.y = TRUE to convert to -log10(y), FALSE to leave as is
# ... is passed onto plot() and company
genome.view.new <- function(x, y, chr, chr.order = as.character(unique(chr)), conv.y = TRUE, ...)
{
  if(conv.y)
    y <- -log10(y)

  #####################
  # Plotting defaults #
  #####################

  args <- list(...)

  if(is.null(args$xaxt))
    args$xaxt <- 'n'

  if(is.null(args$bty))
    args$bty <- 'n'

  if(is.null(args$type))
    args$type <- 'l'

  if(is.null(args$xlab))
    args$xlab <- 'Position in Mb'

  if(is.null(args$ylab))
    args$ylab <- expression(paste(-log[10], " p-value", sep = ''))

  if(is.null(args$ylim))
    args$ylim <- range(y)

  if(is.null(args$col)){
    color <- c('blue', 'black', 'green')[1:dim(y)[2]]
  }else{
    color <- rep(args$col, dim(y)[2])[1:dim(y)[2]]
  }

  if(is.null(args$lty)){
    linetype <- 1:dim(y)[2]
  }else{
    linetype <- rep(args$lty, dim(y)[2])[1:dim(y)[2]]
  }

  if(is.null(args$lwd)){
    linewdth <- rep(1, dim(y)[2])
  }else{
    linewdth <- rep(args$lwd, dim(y)[2])[1:dim(y)[2]]
  }

  ##############################################
  # Get chromosome lengths and genome position #
  ##############################################

  chr.min <- sapply(chr.order, function(i) min(x[chr == i]))
  chr.max <- sapply(chr.order, function(i) max(x[chr == i]))
  chr.len <- chr.max - chr.min

  # add 0 to chr 1, don't add last chromosome's length to anything
  gen.pos.add <- c(0, cumsum(as.numeric(chr.len[-length(chr.len)])))
  names(gen.pos.add) <- as.character(chr.order)

  # make sure this evaluates to a character (not numeric)
  if(is.numeric(chr))
    chr <- as.character(chr)

  # convert chromosome position to genome position
  x.geno <- x - chr.min[chr] + gen.pos.add[chr]

  ###################
  # Plot the figure #
  ###################

  # set up plotting region
  args$x <- NA
  args$xlim <- range(x.geno)
  args$y <- NA
  do.call('plot', args)

  # draw curves
  for(i in 1:dim(y)[2])
  {
    lines(x.geno, y[,i], col = color[i], lty = linetype[i], lwd = linewdth[i])
  }

  # mark SNPs
  axis(side = 1, at = x.geno, tick = TRUE, labels = FALSE)

  # render chromosomes
  ycorrect <- (args$ylim[2] - args$ylim[1])/40
  ymin <- par()$usr[3] + ycorrect
  xrange <- par()$usr[2] - par()$usr[1]

  xcorrect <- c(1, -1)*xrange/250 # to make up for rounded end

  for(i in 1:length(chr.order))
  {
    if(i == 1){
      lines(gen.pos.add[i:(i+1)] + c(0, xcorrect[2]), rep(ymin, 2), lwd=7, col='grey50',
            lend = 'round')
    }else{ if(i == length(chr.order)){
      lines(c(gen.pos.add[i], gen.pos.add[i] + chr.len[i]) + c(xcorrect[1], 0), rep(ymin, 2), lwd=7, col='grey50',
            lend = 'round')
    }else{
      lines(gen.pos.add[i:(i+1)] + xcorrect, rep(ymin, 2), lwd=7, col='grey50',
            lend = 'round')
    }}
  }
}
