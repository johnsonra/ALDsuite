// inference.c
// helper functions for admixture inference of statistical significance
// Randall Johnson
// BSP CCR Genetics Core at Frederick National Laboratory
// SAIC-Frederick, Inc
// Created December 18, 2012
// Last Modified December 18, 2012

#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* #############################################
   # lgs_num(Aj, probs0, probs1, probs2, dims) #
   ############################################# */

// Aj - array of ancestral probabilities
// phi1 - risk associated with 1 risk allele
// phi2 - risk associated with 2 risk alleles
// probs0 - columns of Aj (3rd dim) corresponding to no risk alleles
// probs1 - " " corresponding to 1 risk allele
// probs2 = " " corresponding to 2 risk alleles
// dims = dimensions of Aj
SEXP lgs_num(SEXP Aj, SEXP phi1, SEXP phi2, SEXP probs0, SEXP probs1, SEXP probs2, SEXP dims)
{
  // counters
  int  i,  j,  p,  iter,
      ni, nj, np, niter,
      ngroups; // maximum number of groups of p to sum over

  ni = INTEGER(dims)[0];
  nj = INTEGER(dims)[1];
  np = INTEGER(dims)[2];
  niter = INTEGER(dims)[3];

  ngroups = LENGTH(probs0);
  if(LENGTH(probs1) > ngroups)
  {
    ngroups = LENGTH(probs1);
  }

  // pointers to SEXPs
  double *xAj, *xphi1, *xphi2, *xretval;
  int *xprobs0, *xprobs1, *xprobs2;

  xAj = REAL(Aj);
  xphi1 = REAL(phi1);
  xphi2 = REAL(phi2);
  xprobs0 = INTEGER(probs0);
  xprobs1 = INTEGER(probs1);
  xprobs2 = INTEGER(probs2);

  // sub sum for each marker
  double markerSum;

  // value to be returned
  SEXP retval;
  PROTECT(retval = allocVector(REALSXP, nj*niter));
  xretval = REAL(retval);

  for(j = 0; j < nj*niter; j++)
  {
    xretval[j] = 0; // going to return a sum of the log
  }

  // sum up numerator
  for(iter = 0; iter < niter; iter++)
  {
    for(j = 0; j < nj; j++)
    {
      for(i = 0; i < ni; i++)
      {
	markerSum = 0;

	for(p = 0; p < ngroups; p++)
	{
	  // 0 risk alleles
	  if(LENGTH(probs0) < p)
	  {
	    markerSum = markerSum + xAj[i + j*ni + xprobs0[p]*ni*nj + iter*ni*nj*np];
	  }

	  // 1 risk allele
	  if(LENGTH(probs1) < p)
	  {
	    markerSum = markerSum + xAj[i + j*ni + xprobs1[p]*ni*nj + iter*ni*nj*np] * xphi1[0];
	  }

	  // 2 risk alleles
	  if(p == 0) // there will only ever be one of these
	  {
	    markerSum = markerSum + xAj[i + j*ni + xprobs2[p]*ni*nj + iter*ni*nj*np] * xphi2[0];
	  }
	}

	xretval[j + iter*nj] = xretval[j + iter*nj] + log(markerSum);
      }
    }
  }

  // set dimensions of retval (nj rows, niter columns)
  SEXP dim;
  PROTECT(dim = allocVector(INTSXP, 2));

  INTEGER(dim)[0] = nj;
  INTEGER(dim)[1] = niter;

  setAttrib(retval, R_DimSymbol, dim);

  // clean up
  UNPROTECT(2);

  return retval;
}
