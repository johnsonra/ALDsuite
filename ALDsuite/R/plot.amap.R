# plot.amap.R (ancestrymap)
# Plot an amap object
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# Created May, 2006
# Last Modified April 16, 2008

plot.amap <- function(x, cM = FALSE, add.legend = TRUE, chromosome, ci = 0.95,
                      ci.col = 'grey75', plot.type = 'genome',
                      change.par = TRUE, ...)
{
  #### THIS ONLY PLOTS DISTANCE MEASURED IN Mb ... ANCESTRYMAP DOES cM AS WELL
  # get phys.pos into a number: 0, at beginning of chr 1 to n, at end of chr 23
  chr.len <- c(245.2, 243.3, 199.4, 191.6, 181.0, 170.7, 158.4, 145.9, 134.5,
               135.5, 135.0, 133.5, 114.2, 105.3, 100.1, 90.0, 81.7, 77.8, 63.8,
               63.6, 47.0, 49.5, 152.6)

  # make sure we are dealing with the correct object
  if(class(x) != 'amap')
    stop("x must be an 'amap' object")

  # start an output device if none is active
  if(dev.cur() == 1)
  {
    stop.dev <- TRUE
    pdf(file='Rplots.pdf', height=8, width=10.5, paper='special')
    warning(paste("Saving plot at ", getwd(), "/Rplots.pdf", sep=''))
  }else{
    stop.dev <- FALSE
  }

  # plot requested plots
  if('chromosome' %in% plot.type)
    chromosome.view(x, chr.len, add.legend, ...)

  if('densFun' %in% plot.type)
    densFun.view(x, chromosome, ci, ci.col, ...)

  if('genome' %in% plot.type)
    genome.view(x, cM, add.legend, change.par, chr.len, ...)

  # close device if needed
  if(stop.dev)
    dev.off()
}
