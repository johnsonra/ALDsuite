// helpers.c
// Helper functions for the ALDsuite package
// Randall Johnson
// BSP CCR Genetics Core at Frederick National Laboratory
// SAIC Frederick, Inc
// Created September 10, 2013
// Last Modified November 9, 2013

#include <R.h>
#include <Rinternals.h>
#include "helpers.h"

/*#############################################################################
  # return one rows per column of a matrix, allowing some rows to be repeated #
  #############################################################################*/

SEXP get_rows(SEXP mat, SEXP rows)
{
  double *xmat = REAL(mat);
  int *xrows = INTEGER(rows);
  int ncols = LENGTH(rows);
  int nrows = LENGTH(mat) / ncols;

  SEXP retval;

  PROTECT(retval = allocVector(REALSXP, ncols));

  double *xretval = REAL(retval);

  for(int j=0; j < ncols; j++)
  {
    xretval[j] = xmat[j*nrows + xrows[j]];
  }

  UNPROTECT(1);
  return retval;
}

/* ######################
   # Indexing functions #
   ###################### */

int index2(int  i, int j,
           int ni)
{
  return i + ni*j;
}

int index3(int  i, int  j, int  k,
           int ni, int nj)
{
  return i + ni*j + ni*nj*k;
}

int index4(int  i, int  j, int  k, int  l,
           int ni, int nj, int nk)
{
  return i + ni*j + ni*nj*k + ni*nj*nk*l;
}

int index5(int  i, int  j, int  k, int  l, int m,
           int ni, int nj, int nk, int nl)
{
  return i + ni*j + ni*nj*k + ni*nj*nk*l + ni*nj*nk*nl*m;
}
