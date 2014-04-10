// sample.joint.gammas.c
// C functions used in sample.joint.gammas()
// Randall Johnson
// BSP CCR Genetics core at Frederick National Laboratory for Cancer Research
// Leidos Biomedical Research, Inc
// Created November 9, 2013
// Last Modified November 9, 2013

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "helpers.h"

/* #######
   # P_O #
   ####### */

// O_prev is an array of haplotypes and has 3 dimensions with the second
//        possibly being of size 1
// P is the allele frequency - a scalar
// prior_betas contains the beta vector for the glm
//             (of length ncol of O_prev + 1)
// prior_eig contains the eigenvalue matrix
//           (with rows equal to ncol of O_prev)
// linked is a logical value indicating if there are linked markers we model
// dims is an integer vector with the dims of O_prev
// pO is what will be returned, real matrix with same number of
//    rows as O_prev and 2 columns (one for mom one for dad)
SEXP P_O(SEXP O_prev, SEXP P, SEXP prior_betas, SEXP prior_eig, SEXP linked,
         SEXP dims)
{
  SEXP pO;
  PROTECT(pO = allocVector(REALSXP, INTEGER(dims)[0]*2));

  double *xOprev = REAL(O_prev);
  double *xP = REAL(P);
  double *xbetas = REAL(prior_betas);
  double *xeig = REAL(prior_eig);
  int *xlinked = INTEGER(linked);
  int *xdims = INTEGER(dims);
  double *xpO = REAL(pO);

  int nbetas = LENGTH(prior_betas);
  int nevecs = nbetas - 1;

  // initialize
  for(int i=0; i < xdims[0]*2; i++)
    xpO[i] = xP[0];

  if(xlinked && nevecs > 0) // they'll just get the allele frequency
  {
    for(int i=0; i < xdims[0]; i++)
    {
      for(int c=0; c < 2; c++)
      {
	double pcs[nevecs];

        // haplotype %*% prior$eig
	for(int v=0; v < nevecs; v++)
	{
	  for(int j=0; j < xdims[1]; j++)
	  {
	    double tmp = xOprev[index3(i, j, c, xdims[0], xdims[1])] *
                         xeig[index2(j, v, xdims[1])];
	    if(j == 0)
	      pcs[v] = tmp;
	    else
	      pcs[v] += tmp;
	  }
	}

        // cbind(1, pcs) %*% prior$betas
	double lod;

	for(int b=0; b < nbetas; b++)
	{
	  if(b == 0)
	    lod = xbetas[b];
	  else
	    lod += xbetas[b] * pcs[b-1];
	}

	double odds = exp(lod);

      	if(isfinite(odds))
      	  xpO[index2(i, c, xdims[0])] = odds / (1 + odds);
      	else
      	{
      	  if(isinf(odds))
      	    xpO[index2(i, c, xdims[0])] = 1;
      	  else
      	    xpO[index2(i, c, xdims[0])] = 0;
      	}
      }
    }
  }

  SEXP pO_dims;
  PROTECT(pO_dims = allocVector(INTSXP, 2));

  INTEGER(pO_dims)[0] = xdims[0];
  INTEGER(pO_dims)[1] = 2;

  setAttrib(pO, R_DimSymbol, pO_dims);

  UNPROTECT(2);
  return pO;
}
