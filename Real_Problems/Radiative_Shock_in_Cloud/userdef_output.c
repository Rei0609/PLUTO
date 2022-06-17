#include "pluto.h"
#include "radiat_frac.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k, nv;
  double v[NVAR], rhs[NVAR], em[7];
  double *x1, *x2, *x3;
  double ***te, mu;

  te = GetUserVar("te");

  DOM_LOOP(k, j, i) {

    /* Temperature */
    NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
    mu = MeanMolecularWeight(v);
    te[k][j][i] = v[PRS] / v[RHO] * mu;

  }

#if COOLING == TABULATED

  double ***lmd, ***l_oiii;

  /* New variables - names must exist under uservar */
  lmd = GetUserVar("lmd");
  l_oiii = GetUserVar("l_oiii");

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  DOM_LOOP(k, j, i) {

    NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
    mu = MeanMolecularWeight(v);

    /* Cooling rate */
    v[RHOE] /= g_gamma - 1.;
    Radiat(v, rhs);
    lmd[k][j][i] = rhs[RHOE];

    Radiat_frac(v, em);
    l_oiii[k][j][i] = em[1];
  }

#endif


}

/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{
  Image *image;

  /* HDF5 output cannot be controlled yet. Everything is output.*/

  /* FLT output */
  SetOutputVar("prs",     FLT_OUTPUT, YES);
  SetOutputVar("te",      FLT_OUTPUT, YES);
#if COOLING == TABULATED
  SetOutputVar("lmd",     FLT_OUTPUT, YES);
  SetOutputVar("l_oiii",  FLT_OUTPUT, YES);
#endif


#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}

