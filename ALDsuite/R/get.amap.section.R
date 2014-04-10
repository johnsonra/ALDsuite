# get.amap.section.R (ancestrymap)
# get a section from ANCESTRYMAP output
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created August 2005


get.amap.section <- function(x, section.head, file.out='temp.out.not.used', header=TRUE, messed.up.header=FALSE)
{
  # find start point
  to.skip <- c(1:length(x))[x == section.head][1] # find section header

  if(is.na(to.skip)) # if it isn't found...
    return(NULL)

  # skip following comments
  while((substr(x[to.skip + 1], 1, 1) == '#' & substr(x[to.skip + 1], 1, 3) != '###') | x[to.skip + 1] %in% c("", " "))
    to.skip <- to.skip + 1

  if(substr(x[to.skip + 1], 1, 3) == '###')
    return(NULL) # if there are no data...

  # some headers have one too few headings, lets just ignore the headers
  if(messed.up.header)
    header <- FALSE

  to.stop <- c(1:length(x))[substr(x, 1, 1) == "#"] # find all possible stopping places
  if(!any(to.stop > to.skip))
    to.stop <- length(x)
  else
    to.stop <- min(to.stop[to.stop > to.skip]) # choose the one right after the starting place

  # compensate for blank lines (ignored in read.table statement)
  temp <- x[(to.skip+1):(to.stop-1)]
  n.lines2read <- to.stop - to.skip - ifelse(header, 2, 1) - max(c(length(temp[temp %in% c('', ' ')]), 0))

  # read in dataset and make names nicer
  x <- read.table(file.out, sep='', skip=to.skip, nrows=n.lines2read, fill=TRUE, header=header)
  
  names(x) <- gsub('_', '.', tolower(names(x)))

  if(messed.up.header)
    x <- x[-1,]

  return(x)
}
