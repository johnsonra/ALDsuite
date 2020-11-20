# phys.calc.R
# Calculate physical position given genetic position
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created May 18, 2007
# Last Modified January 4, 2012

phys.calc <- function(chrom, pos, n.extrap = 4)
{
  # check input
  if(length(chrom) != length(pos))
    stop("Chromosome and Genetic positions must be the same length")

  if(any(pos < 1))
    warning('Are your positions in Morgans or Centimorgans?\nShould be cM')
  
  # set warning options back to normal when exiting
  w <- options("warn")
  on.exit(options(w))

  # turn off warnings
  options(warn = -1)
  
  phys.pos <- lapply(1:length(chrom), function(i)
                     try(uniroot(.phys.calc, interval = c(0, 250e6), chrom = chrom[i],
                      pos = pos[i], n.extrap = n.extrap)))
  
  phys.pos <- sapply(phys.pos, function(x){
                     if(class(x) == 'try-error')
                       return(NA)
                     return(x$root)})

  return(data.frame(chr = chrom,
                    phys.pos = phys.pos,
                    gen.pos = pos))
}

.phys.calc <- function(chrom, pos, n.extrap, guess)
{
  gen.calc(chrom, guess, n.extrap)$gen.pos - pos
}
