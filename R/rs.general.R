# rs.general.R
# Get general information about a SNP
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 22, 2006
# Last Modified October 27, 2006

rs.general <- function(xml)
{
  # get the organism -- can't seem to find this in the most recent version :(
##   organism <- trim(grep('&lt;Rs_organism&gt;', xml, value=TRUE))
##   organism <- sub('&lt;/Rs_organism&gt;', '',
##                   sub('&lt;Rs_organism&gt;', '', organism))

  # get the linenumber with the rs id info
  rsnumber <- grep('Rs rsId=', xml)

  # get the sequence information directly after the rs id info
  tag5 <- grep('Seq5', xml)
  tag5 <- trim(xml[min(tag5[tag5 > rsnumber])])
  tag5 <- sub('&lt;Seq5&gt;', '',
          sub('&lt;/Seq5&gt;', '', tag5))
  
  tag3 <- grep('Seq3', xml)
  tag3 <- trim(xml[min(tag3[tag3 > rsnumber])])
  tag3 <- sub('&lt;Seq3&gt;', '',
          sub('&lt;/Seq3&gt;', '', tag3))
  
  polymorph <- grep('Observed', xml)
  polymorph <- trim(xml[min(polymorph[polymorph > rsnumber])])
  polymorph <- sub('&lt;Observed&gt;', '',
               sub('&lt;/Observed&gt;', '', polymorph))
  
  # get line with rs number and drop unwanted bits
  rsnumber <- trim(xml[rsnumber])
  rsnumber <- sub('&lt;Rs ', '',
              sub('&gt;', '',
              gsub('&quot;', '', rsnumber)))

  # separate each component of the line and then make a list of each component
  rsnumber <- strsplit(rsnumber, ' ', fixed = TRUE)[[1]]
  rsnumber <- strsplit(rsnumber, '=', fixed = TRUE)

  # extract names and values of each component
  name <- sapply(rsnumber, function(x) x[1])
  rsnumber <- sapply(rsnumber, function(x) x[2])
  names(rsnumber) <- name

  # get heterozygosity
  het <- trim(grep('Het type=', xml, value=TRUE))
  het <- sub('&lt;Het ', '',
         sub('/&gt;', '',
         gsub('&quot;', '', het)))

  # separate each component of the line and then make a list of each component
  het <- strsplit(het, ' ', fixed = TRUE)[[1]]
  het <- strsplit(het, '=', fixed = TRUE)

  # extract names and values of each component
  name <- sapply(het, function(x) x[1])
  het <- sapply(het, function(x) x[2])
  names(het) <- name

  return(list(#organism = organism,
              rsnumber = paste('rs', rsnumber['rsId'], sep = ''),
              type = rsnumber['snpClass'],
              het = het['value'],
              het.se = het['stdError'],
              tag5,
              polymorph,
              tag3))
}
