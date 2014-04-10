// likelihoods.c
// C functions to calculate likelihoods
// Randall Johnson
// BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
// SAIC-Frederick, Inc
// Created October 31, 2012
// Last Modified December 4, 2013

#include <R.h>
#include <Rinternals.h>
/* #include <RMath.h> */
#include <math.h>


/*############################################
  # P(recombination == x ; lambda, distance) #
  ############################################*/

SEXP P_recombination(SEXP x, SEXP lambda, SEXP distance)
{
  int i, nlambda, ndistance, n, *xx;
  double *xlambda, *xdistance, *xretval;
  SEXP retval;

  // figure out how big n should be
  nlambda = LENGTH(lambda);
  ndistance = LENGTH(distance);

  if(nlambda > ndistance)
    n = nlambda;
  else
    n = ndistance;

  // type checking/casting
  PROTECT(retval = allocVector(REALSXP, n));

  xx = INTEGER(x);
  xlambda = REAL(lambda);
  xdistance = REAL(distance);
  xretval = REAL(retval);

  if(ndistance == nlambda)
  {
    for(i = 0; i < nlambda; i++)
    {
      // P(*no* recombination in lambda generations)
      double noR = pow((1 + exp(-2 * xdistance[i])) / 2, xlambda[i]);

      if(INTEGER(x)[0] == 1)
	xretval[i] =  1 - noR;
      else
	xretval[i] = noR;
    }
  }

  if(ndistance == 1 & n > 1)
  {
    for(i = 0; i < n; i++)
    {
      double noR = pow((1 + exp(-2 * xdistance[0])) / 2, xlambda[i]);

      if(xx[0] == 1)
	xretval[i] = 1 - noR;
      else
	xretval[i] = noR;
    }
  }

  if(nlambda == 1 && n > 1)
  {
    for(i = 0; i < n; i++)
    {
      double noR = pow((1 + exp(-2 * xdistance[i])) / 2, xlambda[0]);

      if(xx[0] == 1)
	xretval[i] = 1 - noR;
      else
	xretval[i] = noR;
    }
  }

  UNPROTECT(1);
  return retval;
}


/*###########################
  # logP(gammas & geno | P) #
  ###########################*/

SEXP logP_gammas_geno_P_j(SEXP gammas1, SEXP gammas2, SEXP geno, SEXP P)
{
  int i, pop1, pop2, nindivs, npops, *xgeno;
  double *xgammas1, *xgammas2, *xP, *xretval;
  SEXP retval;

  xgammas1 = REAL(gammas1);
  xgammas2 = REAL(gammas2);
  xgeno = INTEGER(geno);
  xP = REAL(P);

  npops = LENGTH(P);
  nindivs = LENGTH(geno);

  double ps [nindivs];

  // initialize ps
  for(i = 0; i < nindivs; i++)
  {
    ps[i] = 0;
  }

  // calculate ps
  for(pop1 = 0; pop1 < npops; pop1++)
  {
    for(pop2 = 0; pop2 < npops; pop2++)
    {
      for(i = 0; i < nindivs; i++)
      {
	if(xgeno[i] == 0)
	  ps[i] = ps[i] + (1 - xP[pop1]) * xgammas1[pop1*nindivs + i] *
                          (1 - xP[pop2]) * xgammas2[pop2*nindivs + i];

	if(xgeno[i] == 1)
	  ps[i] = ps[i] + (1 - xP[pop1]) * xgammas1[pop1*nindivs + i] *
                               xP[pop2]  * xgammas2[pop2*nindivs + i] / 2 +
                               xP[pop1]  * xgammas1[pop1*nindivs + i] *
                          (1 - xP[pop2]) * xgammas2[pop2*nindivs + i] / 2;

	if(xgeno[i] == 2)
	  ps[i] = ps[i] + xP[pop1] * xgammas1[pop1*nindivs + i] *
	                  xP[pop2] * xgammas2[pop2*nindivs + i];
      }
    }
  }

  PROTECT(retval = allocVector(REALSXP, 1));
  xretval = REAL(retval);

  // sum up ps
  for(i = 0; i < nindivs; i++)
  {
    if(i == 0)
      xretval[0] = log(ps[0]);
    else
      xretval[0] = xretval[0] + log(ps[i]);
  }

  UNPROTECT(1);

  return retval;
}


/*###########################
  # logP(G & geno | P) #
  ###########################*/

SEXP logP_Gs_geno_P_j(SEXP G1, SEXP G2, SEXP geno, SEXP P)
{
  int i, l, npops, nindivs,
      *xG1, *xG2, *xgeno;
  SEXP retval;
  double ps, *xP, *xretval; // ps will store the individual probabilities

  // P for each individual based on G, and Q = (1 - P) as well
  double P1, P2, Q1, Q2;

  xG1 = INTEGER(G1);
  xG2 = INTEGER(G2);
  xgeno = INTEGER(geno);
  xP = REAL(P);

  nindivs = LENGTH(geno);
  npops = LENGTH(P);

  PROTECT(retval = allocVector(REALSXP, 1));
  xretval = REAL(retval);
  xretval[0] = 0;

  // do calculations on a per individual bases and keep a running talley in retval :)
  for(i = 0; i < nindivs; i++)
  {
    // initialize each value for each individual
    P1 = 0;
    P2 = 0;
    Q1 = 0;
    Q2 = 0;

    for(l = 0; l < npops; l++)
    {
      P1 = P1 + xG1[i + l*nindivs] * xP[l];
      P2 = P2 + xG2[i + l*nindivs] * xP[l];
      Q1 = Q1 + xG1[i + l*nindivs] * (1 - xP[l]);
      Q2 = Q2 + xG2[i + l*nindivs] * (1 - xP[l]);
    }

    if(xgeno[i] == 0)
      ps = Q1*Q2;

    if(xgeno[i] == 1)
      ps = Q1*P2 + P1*Q2;

    if(xgeno[i] == 2)
      ps = P1*P1;

    xretval[0] = xretval[0] + log(ps);
}

  UNPROTECT(1);

  return retval;
}

/*####################################
  # P(a == x | gamma, gamma.alt ; P) #
  ####################################*/


SEXP P_a_gamma_gamma1(SEXP x, SEXP gamma, SEXP gamma_alt, SEXP P)
{
  int i, n, *xx, *xgamma, *xgamma_alt;
  double *xP, *xretval;

  SEXP retval;

  xx = INTEGER(x);
  xgamma = INTEGER(gamma);
  xgamma_alt = INTEGER(gamma_alt);
  xP = REAL(P);

  n = LENGTH(x);

  PROTECT(retval = allocVector(REALSXP, n));
  xretval = REAL(retval);

  for(i = 0; i < n; i++)
  {
    switch(xx[i]){
    case 0:
      xretval[i] = (1 - xP[xgamma[0]]) * (1 - xP[xgamma_alt[0]]);
      break;

    case 1:
      xretval[i] = xP[xgamma[0]] * (1 - xP[xgamma_alt[0]]) +
	           (1 - xP[xgamma[0]]) * xP[xgamma_alt[0]];
      break;

    case 2:
      xretval[i] = xP[xgamma[0]] * xP[xgamma_alt[0]];
      break;

    default: // if missing return 1
      xretval[i] = 1;
      break;
    }
  }

  UNPROTECT(1);

  return(retval);
}

SEXP P_a_gamma_gamma2(SEXP x, SEXP gamma, SEXP gamma_alt, SEXP P, SEXP dims)
{
  int i, n, *xx, *xgamma, *xgamma_alt, *xdims;
  double *xP, *xretval;

  SEXP retval;

  xx = INTEGER(x);
  xgamma = INTEGER(gamma);
  xgamma_alt = INTEGER(gamma_alt);
  xP = REAL(P);
  xdims = INTEGER(dims);

  int pop = xdims[0] * xdims[1]; // offset for population (3rd dimension of P)
  int father = xdims[0]; // offset for father's chromosome (2nd dimension of P)

  n = LENGTH(x);

  PROTECT(retval = allocVector(REALSXP, n));
  xretval = REAL(retval);

  for(i = 0; i < n; i++)
  {
    switch(xx[i]){
    case 0:
      xretval[i] = (1 - xP[pop*xgamma[0] + i]) * (1 - xP[pop*xgamma_alt[0] + father + i]);
      break;

    case 1:
      xretval[i] =      xP[pop*xgamma[0] + i] * (1 - xP[pop*xgamma_alt[0] + father + i]) +
	           (1 - xP[pop*xgamma[0] + i]) *     xP[pop*xgamma_alt[0] + father + i];
      break;

    case 2:
      xretval[i] = xP[pop*xgamma[0] + i] * xP[pop*xgamma_alt[0] + father + i];
      break;

    default: // if missing return 1
      xretval[i] = 1;
      break;
    }
  }

  UNPROTECT(1);

  return(retval);
}


/* ###########################
   # logP(beta | sigma, tau) #
   ########################### */

// to be called from C
double logP_beta_sigma_tau(int n, double* xbeta, double* xbeta_p, double* xsigma, double* xtau)
{
  double retval = n*log(xtau[0]);

  // retval - sum(log(sigma) + [tau * (beta - beta_p)]^2 / sigma)
  for(int i=0; i < n; i++)
  {
    retval -= log(xsigma[i]) + pow(xtau[0] * (xbeta[i] - xbeta_p[i]), 2) / xsigma[i];
  }

  retval = retval / 2;

  return retval;
}

// to be called from R
SEXP logP_beta_sigma_tauR(SEXP beta, SEXP beta_p, SEXP sigma, SEXP tau)
{
  SEXP retval;
  PROTECT(retval = allocVector(REALSXP, 1));

  REAL(retval)[0] = logP_beta_sigma_tau(LENGTH(beta), REAL(beta), REAL(beta_p), REAL(sigma), REAL(tau));

  UNPROTECT(1);
  return retval;
}


/* SEXP sample_tau_Pm_prior(SEXP beta, SEXP beta_p, SEXP beta_prior, SEXP sigma, SEXP tau, SEXP tau_new) */
/* { */
/*   // number of parameters / dimensions of n x n square matrix, hessian */
/*   int n = LENGTH(beta); */

/*   // retval */
/*   SEXP retval; */
/*   PROTECT(retval = allocVector(REALSXP, n + 2)); // two logliks + one for each beta */

/*   // pointers to SEXPs */
/*   double *xbeta = REAL(beta); */
/*   double *xbeta_p = REAL(beta_p); */
/*   double *xbeta_prior = REAL(beta_prior); */
/*   double *xsigma = REAL(sigma); */
/*   double *xtau = REAL(tau); */
/*   double *xtau_new = REAL(tau_new); */
/*   double *xretval = REAL(retval); */

/*   // get likelihood of current configuration -- if I store this, I could save some time? */
/*   double loglik_cur = logP_beta_sigma_tau(n, xbeta, xbeta_p, xsigma, xtau); */

/*   // new beta_p */
/*   double new_beta_p [n]; */

/*   GetRNGstate(); */
/*   for(int l=0; l < n; l++) */
/*   { */
/*     new_beta_p[l] = rnorm(xbeta_prior[l], xsigma[l]); */
/*   } */

/*   // likelihood of new beta */
/*   double loglik_new = logP_beta_sigma_tau(n, xbeta, (&new_beta_p)[n], xsigma, xtau); */

/*   // decision */
/*   double p_switch = exp(loglik_new - loglik_cur); */

/*   double decision_point = rand(); */
/*   PutRNGstate(); */

/*   // save likelihoods and updated value for beta */
/*   if(p_switch > decision_point) */
/*   { */
/*     xretval[0] = loglik_new; */
/*     xretval[1] = logP_beta_sigma_tau(n, xbeta, (&new_beta_p)[n], xsigma, xtau_new); */

/*     for(int l=0; l < n; l++) */
/*     { */
/*       xretval[l + 2] = new_beta_p[l]; */
/*     } */
/*   }else{ */
/*     xretval[0] = loglik_cur; */
/*     xretval[1] = logP_beta_sigma_tau(n, xbeta, xbeta_p, xsigma, xtau_new); */

/*     for(int l=0; l < n; l++) */
/*     { */
/*       xretval[l + 2] = xbeta_p[l]; */
/*     } */
/*   } */

/*   UNPROTECT(1); */

/*   return retval; */
/* } */
