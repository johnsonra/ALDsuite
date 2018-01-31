# hCV.lookup.R
# Get position and chromosome of a SNP (given a Celera number)
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 21, 2006
# Last Modified October 27, 2008

hCV.lookup <- function(hCVnumber)
{
  # make sure we are dealing with the correct data type
  if(!is.character(hCVnumber))
    hCVnumber <- as.character(hCVnumber)

  # check to be sure we have a good number
  if(substr(hCVnumber, 1, 3) != 'hCV' |
     is.na(as.numeric(sub('hCV', '', hCVnumber))))
  {
    warning(paste(hCVnumber, "seems to be an invalid hCV number"))
    return(NULL)
  }

  # get html code from NCBI
  html <- try(readLines(paste('http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Snp&term=',
                              hCVnumber, sep=''), warn = FALSE), silent = TRUE)

  # if internet connection is down (either locally or at NCBI)
  if(class(html) == 'try-error')
    return(NULL)
    
  # get the celera number if one exists
  html <- grep('snp_ref.cgi', html, value=TRUE)

  if(length(html) == 0)
    return(NULL)

  html <- unlist(strsplit(html, '\">rs'))[2]
  html <- unlist(strsplit(html, '</a>'))[1]

  # get the rs number information
  x <- rs.lookup(paste('rs', html, sep=''))
  x$hCVnumber <- hCVnumber

  return(x)
}
