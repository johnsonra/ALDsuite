# dbSNP.R
# Look up information in dbSNP
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 22, 2006
# Last Modified October 27, 2008

dbSNP <- function(snp, parallel = FALSE, processes = 2, type = 'MPI')
{
##   if(parallel)
##   {
##     if(!require(snow))
##       stop('snow required for the parallel option... see ?dbSNP')
    
##     cl <- makeCluster(processes, type=type)

##     clusterEvalQ(cl, library(gdata))

##     snps <- clusterApplyLB(cl, snp, lookup)
    
##     stopCluster(cl)
##   }else{
    snps <- lapply(snp, lookup)
##   }

  snps <- data.frame(organism = get.snp.part(snps, 'organism'),
                     rsnumber = get.snp.part(snps, 'rsnumber'),
                     type = get.snp.part(snps, 'type'),
                     het = get.snp.part(snps, 'het'),
                     het.se = get.snp.part(snps, 'het.se'),
                     tag5 = get.snp.part(snps, 'tag5'),
                     polymorph = get.snp.part(snps, 'polymorph'),
                     tag3 = get.snp.part(snps, 'tag3'),
                     genomeBuild = get.snp.part(snps, 'genomeBuild'),
                     dbSnpBuild = get.snp.part(snps, 'dbSnpBuild'),
                     chr = get.snp.part(snps, 'chr'),
                     phys.pos = get.snp.part(snps, 'phys.pos') + 1, # I'm not sure why NCBI is 1 off, but it is...
                     hCVnumber = get.snp.part(snps, 'hCVnumber'),
                     stringsAsFactors = FALSE
                     )

  return(snps)
}

lookup <- function(snp)
{
  if(substr(snp, 1, 1) == 'r'){
    return(rs.lookup(snp))
  }else{ if(substr(snp, 1, 1) == 'h'){
    return(hCV.lookup(snp))
  }else{
    warning(paste('Bad SNP ', snp, '?', sep=''))
    return(NULL)
  }}
}

  
get.snp.part <- function(snps, part)
{
  unlist(sapply(snps, function(x)
                {
                  if( is.null(x[[part]]) | length(x[[part]]) == 0){
                    return(NA)
                  }else{
                    return(x[[part]])
                  }
                }
               )
        )
}

