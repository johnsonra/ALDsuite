# rs.lookup.R
# Get position and chromosome of a SNP
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 20, 2006
# Last Modified October 27, 2008

rs.lookup <- function(rsnumber)
{
  # make sure we are dealing with the correct data type
  if(!is.character(rsnumber))
    rsnumber <- as.character(rsnumber)

  # check to be sure we have a good number
  if(substr(rsnumber, 1, 2) != 'rs' | is.na(as.numeric(sub('rs', '', rsnumber))))
  {
    warning(paste(rsnumber, 'seems to be an invalid rs number'))
    return(NULL)
  }

  # get xml code from NCBI
  xml <- try(readLines(paste('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&report=XML&id=',
                             sub('rs', '', rsnumber), sep=''), warn = FALSE))

  # if internet connection is down (either locally or at NCBI)
  if(class(xml) == 'try-error')
    return(NULL)
  
  # Warn and return NULL if the SNP is not found at NCBI
  if(xml[1] == "<ERROR>Empty id list - nothing todo</ERROR>" |
     as.logical(length(grep('Error occurred', xml))))
  {
    warning(paste(rsnumber, "NCBI error or Not found at NCBI\n"))
    return(NULL)
  }

  # get general snp info
  general <- rs.general(xml)
  
  # get information from Assembly chunk
  assembly <- rs.assembly(xml)

  # note that there may be other components in the current reference assembly,
  # there should only be one where reference == TRUE
  x <- list(#organism = general$organism, # find this! :( :(
            rsnumber = general$rsnumber,
            type = general$type,
            het = as.numeric(general$het),
            het.se = as.numeric(general$het.se),
            tag5 = general$tag5,
            polymorph = general$polymorph,
            tag3 = general$tag3,
            genomeBuild = assembly$genomeBuild,
            dbSnpBuild = assembly$dbSnpBuild,
            chr = assembly$chr,
            phys.pos = as.numeric(assembly$phys.pos))

  return(x)
}
