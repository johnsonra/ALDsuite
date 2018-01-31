// helpers.h
// Function headers for helper functions
// Randall Johnson
// BSP CCR Genetics Core at Frederick National Laboratory
// Leidos Biomedical Research, Inc
// Created November 9, 2013
// Last Modified November 9, 2013


// indexing function for SEXPs
int index2(int  i, int j,
           int ni);

int index3(int  i, int  j, int  k,
           int ni, int nj);

int index4(int  i, int  j, int  k, int  l,
           int ni, int nj, int nk);

int index5(int  i, int  j, int  k, int  l, int m,
           int ni, int nj, int nk, int nl);

// returns one row per column of a matrix, allowing some rows to be repeated
SEXP get_rows(SEXP mat, SEXP rows);
