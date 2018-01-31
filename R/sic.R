# sic.R
# Calcuation of Shannon Information Content
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created December 7, 2006
# Last Modified December 8, 2006

sic <- function(fa, na, fe, ne, admix = .2, type='sic')
{
  # correct for any 0's
  fa <- ifelse(fa == 0, 1 / (na + 1), fa)
  fe <- ifelse(fe == 0, 1 / (ne + 1), fe)
  
  # expected african contribution
  p11 <- (1-admix) * fa
  p21 <- (1-admix) * (1-fa)

  # expected european contribution
  p12 <- admix * fe
  p22 <- admix * (1-fe)

  # expected allele frequencies
  a1 <- p11 + p12
  a2 <- p21 + p22

  ##### calculations ####
  # ∑_i∑_j p_{ij} * log(p_{ij})
  sisj.plp <- p11 * log(p11) + p12 * log(p12) +
              p21 * log(p21) + p22 * log(p22)

  # ∑_i p_i * log(p_i)
  si.plp <- admix * log(admix) + (1-admix) * log(1-admix)

  # ∑_j p_j * log(p_j)
  sj.plp <- a1 * log(a1) + a2 * log(a2)

  #∑_i∑_j ( p_{ij}^2 * log(p_{ij}) ) /
  #            p_{ij} - p_{i'j}
  sisj.p2lp.dpi <- (p11^2 * log(p11)) / (p11 - p21) +
                   (p12^2 * log(p12)) / (p12 - p22) +
                   (p21^2 * log(p21)) / (p21 - p11) +
                   (p22^2 * log(p22)) / (p22 - p12)
  
  # Information
  In <- switch(type,
                 sic = sisj.plp - si.plp - sj.plp,
                 In = sisj.plp / 2 - sj.plp,
                 Ia = sisj.p2lp.dpi - 1 - sj.plp,
                 NA)

  return(In)
}
