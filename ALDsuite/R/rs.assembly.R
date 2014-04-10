# rs.assembly.R
# Get information from the Assembly chunk
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 22, 2006
# Last Modified October 27, 2008

rs.assembly <- function(xml)
{
  # get the starting line of the current reference assembly build
  assembly <- grep('reference=&quot;true&quot;', xml)

  # now, given the start line of the correct assembly chunk,
  # get chrmosome number
  chr <- grep('Component componentType', xml)
  chr <- trim(xml[min(chr[chr > assembly])])
  chr <- sub('&lt;Component ', '',
         sub('&gt;', '',
         gsub('&quot;', '', chr)))

  # separate each component
  chr <- strsplit(chr, ' ', fixed = TRUE)[[1]]
  chr <- strsplit(chr, '=', fixed = TRUE)

  # extract names and values
  name <- sapply(chr, function(x) x[1])
  chr <- sapply(chr, function(x) x[2])
  names(chr) <- name

  # now, given the start line of the correct assembly chunk,
  # get physical position
  phys.pos <- grep('MapLoc', xml)
  phys.pos <- trim(xml[min(phys.pos[phys.pos > assembly])])
  phys.pos <- sub('&lt;MapLoc ', '',
              sub('&gt;&lt;/MapLoc&gt;', '',
              gsub('&quot;', '', phys.pos)))

  # separate each component
  phys.pos <- strsplit(phys.pos, ' ', fixed = TRUE)[[1]]
  phys.pos <- strsplit(phys.pos, '=', fixed = TRUE)
  
  # extract names and values
  name <- sapply(phys.pos, function(x) x[1])
  phys.pos <- sapply(phys.pos, function(x) x[2])
  names(phys.pos) <- name

  # now, finally, get the assembly information :)
  assembly <- trim(xml[assembly])
  assembly <- sub('&lt;Assembly ', '',
              sub('&gt;', '',
              gsub('&quot;', '', assembly)))

  # separate each component
  assembly <- strsplit(assembly, ' ', fixed = TRUE)[[1]]
  assembly <- strsplit(assembly, '=', fixed = TRUE)

  # extract names and values
  name <- sapply(assembly, function(x) x[1])
  assembly <- sapply(assembly, function(x) x[2])
  names(assembly) <- name
  
  return(list(genomeBuild = assembly['genomeBuild'],
              dbSnpBuild = assembly['dbSnpBuild'],
              chr = chr['chromosome'],
              phys.pos = phys.pos['physMapInt']))              
}
