// samples.c
// C functions used in samples.R
// Randall Johnson
// BSP CCR Genetics core at Frederick National Laboratory for Cancer Research
// SAIC-Frederick, Inc
// Created December 13, 2012
// Last Modified September 16, 2013

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "helpers.h"

/*###########################
  # count.Os(j, geno, P, G) #
  ###########################*/

SEXP count_Os(SEXP geno, SEXP P, SEXP G1, SEXP G2)
{
  // counters for individual, parent, population and ref/var allele - also total number of individuals and populations
  int i, j, k, l, m, ni, nl;

  ni = LENGTH(geno);
  nl = LENGTH(P);

  // population picked for each parent
  int pop [2];

  // variables used to decide which parent to allocate alleles to
  double decision_point;

  // R objects
  int *xgeno, *xG1, *xG2, *xcounts;
  double *xP;

  xgeno = INTEGER(geno);
  xP = REAL(P);
  xG1 = INTEGER(G1);
  xG2 = INTEGER(G2);

  // initialize value to be returned
  SEXP counts, dim;

  PROTECT(counts = allocVector(INTSXP, ni*nl*2)); // individual, population, ref/vart
  xcounts = INTEGER(counts);

  for(i = 0; i < ni*nl*2; i++)
  {
    xcounts[i] = 0;
  }

  // set up random number generator
  GetRNGstate();

  for(i = 0; i < ni; i++)
  {
    // identify which population was picked for each parent
    for(l = 0; l < nl; l++)
    {
      pop[0] = xG1[ni*l + i] * l; // using 0-based indices on purpose!
      pop[1] = xG2[ni*l + i] * l;
    }

    // assign allele counts (ref/vart) to populations
    switch(xgeno[i]){
    case 0: // if 0, then both are variant alleles
      xcounts[nl*ni*1 + ni*pop[0] + i] = 1;
      xcounts[nl*ni*1 + ni*pop[1] + i] = xcounts[nl*ni*1 + ni*pop[1] + i] + 1;
      break;

    case 1: // if 1, then we have to decide
      decision_point = xP[pop[0]] * (1 - xP[pop[1]]) /
	xP[pop[0]] * (1 - xP[pop[1]]) + xP[pop[1]] * (1 - xP[pop[0]]);

      if(unif_rand() < decision_point){
        // less than decision_point -> parent0 gets the reference allele
	xcounts[nl*ni*0 + ni*pop[0] + i] = 1;
	xcounts[nl*ni*1 + ni*pop[1] + i] = xcounts[nl*ni*1 + ni*pop[1] + i] + 1;
      }else{
        // grater than decision_point -> parent1 get the reference allele
	xcounts[nl*ni*0 + ni*pop[1] + i] = 1;
	xcounts[nl*ni*1 + ni*pop[0] + i] = xcounts[nl*ni*1 + ni*pop[0] + i] + 1;
      }
      break;

    case 2: // if 2, then we have two reference alleles
      xcounts[nl*ni*0 + ni*pop[0] + i] = 1;
      xcounts[nl*ni*0 + ni*pop[1] + i] = xcounts[nl*ni*0 + ni*pop[1] + i] + 1;
      break;

    default: // if missing or something else unexpected, return 0 (already set)
      break;
    }
  }

  // set array dimensions properly before returning
  // catch dimnames in the R function -- too messy here (see pg 111 of R-exts)
  PROTECT(dim = allocVector(INTSXP, 3));

  INTEGER(dim)[0] = ni;
  INTEGER(dim)[1] = nl;
  INTEGER(dim)[2] = 2;

  setAttrib(counts, R_DimSymbol, dim);

  // clean up
  PutRNGstate();
  UNPROTECT(2);

  return counts;
}

/* ##########################################################
   # count_Os(int dims, int geno, double P, int G1, int G2) #
   ########################################################## */

// index in {0,1,2,...}
// number of indicies for the respective index in {1,2,3,...}
//    i = individual
//    j = marker
//    k = parent
//    l = population
//    m = reference/varinat

SEXP sample_O(SEXP dims, SEXP geno, SEXP P, SEXP G)
{
  // R objects
  int *xdims, *xgeno, *xG, *xcounts;
  double *xP;

  xdims = INTEGER(dims);
  xgeno = INTEGER(geno);
  xP = REAL(P);
  xG = INTEGER(G);

  // counters for individual, parent, population and ref/var allele - also total number of individuals and populations
  int i, j, k, l, m, ni, nj, nk, nl;

  ni = xdims[0];
  nj = xdims[1];
  nk = xdims[2]; // this one should always be 2, but include for completeness
  nl = xdims[3];

  // population picked for each parent and index of P for those populations
  int pop[2], p0, p1;

  // variables used to decide which parent to allocate alleles to
  double decision_point;

  // initialize value to be returned
  SEXP counts, newdim;

  PROTECT(counts = allocVector(INTSXP, ni*nj*nl*2)); // individual, marker, population, ref/vart
  xcounts = INTEGER(counts);

  for(i = 0; i < ni*nj*nl*2; i++)
  {
    xcounts[i] = 0;
  }

  // set up random number generator
  GetRNGstate();

  for(i = 0; i < ni; i++)
  {
    for(j = 0; j < nj; j++)
    {
      // identify which population was picked for each parent
      for(l = 0; l < nl; l++)
      {
	pop[0] = xG[index5(i, j, 0, l, 0, ni, nj, nk, nl)] * l; // using 0-based indices on purpose!
	pop[1] = xG[index5(i, j, 1, l, 0, ni, nj, nk, nl)] * l;
      }

      // assign allele counts (ref/vart) to populations
      switch(xgeno[i + ni*j]){
      case 0: // if 0, then both are variant alleles (m = 1)
	xcounts[index5(i, j, 0, pop[0], 1, ni, nj, 1, nl)] = 1;
	m = index5(i, j, 0, pop[1], 1, ni, nj, 1, nl);
	xcounts[m] += 1; // sum for when both parents have same population
	break;

      case 1: // if 1, then we have to decide
	// start by setting these indicies
	p0 = j + pop[0]*nj;
	p1 = j + pop[1]*nj;

	decision_point = xP[p0] * (1 - xP[p1]) /
	  xP[p0] * (1 - xP[p1]) + xP[p1] * (1 - xP[p0]);

	if(unif_rand() < decision_point){
	  // less than decision_point -> parent0 gets the reference allele
	  xcounts[index5(i, j, 0, pop[0], 0, ni, nj, 1, nl)] = 1;
	  m = index5(i, j, 0, pop[1], 1, ni, nj, 1, nl);
	  xcounts[m] = xcounts[m] + 1;
	}else{
          // grater than decision_point -> parent1 get the reference allele
	  xcounts[index5(i, j, 0, pop[1], 0, ni, nj, 1, nl)] = 1;
	  m = index5(i, j, 0, pop[0], 1, ni, nj, 1, nl);
	  xcounts[m] = xcounts[m] + 1;
	}
	break;

      case 2: // if 2, then we have two reference alleles
	xcounts[index5(i, j, 0, pop[0], 0, ni, nj, 1, nl)] = 1;
	m = index5(i, j, 0, pop[1], 0, ni, nj, 1, nl);
	xcounts[m] = xcounts[m] + 1;
	break;

      default: // if missing or something else unexpected, return 0 (already set)
	break;
      }
    }
  }

  // set array dimensions properly before returning
  // catch dimnames in the R function -- too messy here (see pg 111 of R-exts)
  PROTECT(newdim = allocVector(INTSXP, 4));

  INTEGER(newdim)[0] = ni;
  INTEGER(newdim)[1] = nj;
  INTEGER(newdim)[2] = nl;
  INTEGER(newdim)[3] = 2;

  setAttrib(counts, R_DimSymbol, newdim);

  // clean up
  PutRNGstate();
  UNPROTECT(2);

  return counts;
}

/*#####################
  # sample_crossovers #
  #####################*/


SEXP sample_crossovers(SEXP G, SEXP A, SEXP lambda, SEXP d, SEXP chr)
{
  // counters and totals
  int i, j, l, ni, nj, nl;

  ni = LENGTH(lambda);
  nj = LENGTH(d);
  nl = LENGTH(A) / ni;

  // temporary values to be used
  double P0crossovers_marg, Precomb, A1, P0crossovers, q_crossovers, Ptmp, lfact, ld, lld;
  int Gcurr, Gprev, quant, x, odd, chrcurr, chrprev;

  Gcurr = 0;
  chrprev = 0;
  Precomb = 0;

  // this will be the number of crossovers for each person
  SEXP n_crossovers;
  PROTECT(n_crossovers = allocVector(INTSXP, ni));

  // pointers to R objects
  int *xG, *xchr, *xcrossovers;
  double *xA, *xlambda, *xd;

  xcrossovers = INTEGER(n_crossovers);
  xG = INTEGER(G);
  xA = REAL(A);
  xlambda = REAL(lambda);
  xd = REAL(d);
  xchr = INTEGER(chr);

  for(i = 0; i < ni; i++)
  {
    xcrossovers[i] = 0;
  }

  // set up random number generator
  GetRNGstate();

  for(i = 0; i < ni; i++)
  {
    // initialize the number of crossovers
    xcrossovers[i] = 0;


    for(j = 0; j < nj; j++)
    {
      // set chrcurr and chrprev
      chrprev = chrcurr;
      chrcurr = xchr[i];

      // set Gprev to last value of Gcurr
      Gprev = Gcurr;
      Gcurr = 0;

      for(l = 0; l < nl; l++)
      {
      	Gcurr += xG[ni*nj*l + ni*j + i] * l; // population number in [1:nl]
      }

      // skip the rest on the first marker of the chromosome
      if(chrprev == chrcurr & j > 0)
      {
      	// key probabilities
      	P0crossovers_marg = exp(-xlambda[i] * xd[j]);
      	Precomb = (1 - exp(-2 * xlambda[i] * xd[j])) / 2;

      	A1 = xA[ni*Gcurr + i]; // global ancestry for Gcurr population

      	P0crossovers = P0crossovers_marg / ((1 - Precomb) + Precomb*A1);

      	// decide if we observed 1 or more crossovers and count them up
      	if(unif_rand() > P0crossovers || Gcurr != Gprev)
      	{
      	  // decide if there was a recombination event (observed or unobserved)
      	  if(unif_rand() <= Precomb || Gcurr != Gprev)
      	  {
      	    q_crossovers = unif_rand() * Precomb; // q in (0, Precomb)
      	    odd = 1; // odd number of crossovers
      	  }else{
      	    q_crossovers = unif_rand() * (1 - Precomb - P0crossovers_marg);
      	    odd = 0; // even number of crossovers
      	  }

      	  // here is where we count up crosovers
      	  quant = 1 - odd; // start at 0 if odd, 1 if even
      	  Ptmp = 0;
	  lfact = 0;
	  ld = xlambda[i] * xd[j];
	  lld = log(ld);

      	  while(Ptmp < q_crossovers & quant < 10000)
      	  {
      	    x = quant * 2 + odd;

	    // add next two values to lfact (our log x! variable)
	    // we add two, because we are doing every-other value of x
	    lfact += log(x);
	    if(x > 1)
	      lfact += log(x - 1);

      	    Ptmp += exp(x*lld - ld - lfact); // poisson probability on log scale
      	    quant++;
      	  }

      	  xcrossovers[i] += x;
	}
      }
    }
  }

  // clean up
  PutRNGstate();
  UNPROTECT(1);

  return n_crossovers;
}


void sample_G_j(int *G, int j, double *xgammas_joint, int *xchr, int *xmale,
                int *xsex_chr, int *xdims)
{
  int ncombos = xdims[2] * xdims[2];
  int g, k1[ncombos], k2[ncombos];

  // initialize values for population of mother (k1) and father (k2)
  for(int l=0; l < xdims[2]; l++)
  {
    for(int i=0; i < xdims[2]; i++)
    {
      k1[i + l*xdims[2]] = l;
      k2[i + l*xdims[2]] = i;
    }
  }

  GetRNGstate();

  for(int i=0; i < xdims[0]; i ++)
  {
    g = 0;
    double decision_point = unif_rand();
    double p = 1 - xgammas_joint[i + j*xdims[0]]; // g = 0

    // should never need this second catch as sum(gammas.joint[i,j,]) == 1
    // if it isn't, this code will exit with an error rather than entering an infinite loop
    while(decision_point < p && g < ncombos - 1)
    {
      g++;
      p -= xgammas_joint[i + j*xdims[0] + g*xdims[1]*xdims[0]];
    }

    G[i] = k1[g];

    if(xchr[j] == xsex_chr[0] && xmale[i])
    {
      G[i + xdims[0]] = -1;
    }else{
      G[i + xdims[0]] = k2[g];
    }
  }

  PutRNGstate();
}


SEXP sample_G(SEXP gammas_joint, SEXP O, SEXP chr, SEXP male, SEXP sex_chr, SEXP dims)
{
  double *xgammas_joint;
  int *xO, *xchr, *xmale, *xsex_chr, *xdims;

  xgammas_joint = REAL(gammas_joint);
  xO = INTEGER(O);
  xchr = INTEGER(chr);
  xmale = INTEGER(male);
  xsex_chr = INTEGER(sex_chr);
  xdims = INTEGER(dims);

  int G[xdims[0]][2]; // indexed by individual and parent

  SEXP retval;
  // we are returning two arrays here, one for sampled G and one for sampled O
  PROTECT(retval = allocVector(INTSXP, xdims[0]*xdims[1]*2*xdims[2]*2)); // individual, markers, population,
                                                                         // parent/ref/vart, O/G
  int n = LENGTH(retval);
  int *xretval = INTEGER(retval);

  //initiate retval
  for(int i=0; i < n; i++)
  {
    xretval[i] = 0;
  }

  // offsets for different dimensions of retval
  int j_offset = xdims[0]; // marker
  int c_offset = j_offset * xdims[1]; // chromosome/parent or ref/vart allele
  int k_offset = c_offset * 2; // population
  int O_offset = k_offset * xdims[2]; // array of sampled alleles (vs array of sampled ancestral state)

  for(int j=0; j < xdims[1]; j++)
  {
    // reinitialize G before running sample_G_j
    for(int i=0; i < xdims[0]; i++)
    {
      G[i][0] = 0;
      G[i][1] = 0;
    }

    // pick the ancestral state for each individual at each marker
    sample_G_j((int *)G, j, xgammas_joint, xchr, xmale, xsex_chr, xdims);

    for(int i=0; i < xdims[0]; i++)
    {
      int k1 = G[i][0];
      int k2 = G[i][1];
      int g = k1*xdims[2] + k2; // g11, g12, ...

      xretval[i + j*j_offset +            k1*k_offset] = i;//1;
      xretval[i + j*j_offset + c_offset + k2*k_offset] = i;//1;

      for(int k=0; k < xdims[2]; k++)
      {
      	int xIndex = i + j*j_offset + c_offset + k*k_offset + O_offset;
      	int OIndex = i + j*xdims[0] + k*xdims[0]*xdims[1] + g*xdims[0]*xdims[1]*xdims[2];

      	/* xretval[xIndex] = xO[OIndex]; // # variant allele */

      	// count number of chromosomes we sampled from population k
      	int nObs = 0;
      	if(k == k1)
      	  nObs++;
      	if(k == k2)
      	  nObs++;

      	// number of reference alleles is total number of chromosomes minus number of variant alleles
	// remember we are working with 2*O to keep 0.5's integers
      	/* xretval[xIndex - c_offset] = nObs*2 - xO[OIndex]; */
      }
    }
  }

  SEXP newdim;
  PROTECT(newdim = allocVector(INTSXP, 5));

  INTEGER(newdim)[0] = xdims[0];
  INTEGER(newdim)[1] = xdims[1];
  INTEGER(newdim)[2] = 2;
  INTEGER(newdim)[3] = xdims[2];
  INTEGER(newdim)[4] = 2;

  setAttrib(retval, R_DimSymbol, newdim);

  // clean up
  UNPROTECT(2);

  return retval;
}
