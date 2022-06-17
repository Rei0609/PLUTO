/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "init_tools.h"
#include "shock.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  static int once0 = 0;

  /* Some initializations */
  if (!once0) {
    SetBaseNormalization();
    SetIniNormalization();
    once0 = 1;
  }

  /* Input parameters */
  double rho1 = g_inputParam[PAR_RHO1] * ini_code[PAR_RHO1];
  double rho2 = g_inputParam[PAR_RHO2] * ini_code[PAR_RHO2];
  double te2 = g_inputParam[PAR_TE2] * ini_code[PAR_TE2];
  double p1p2 = g_inputParam[PAR_P1P2] * ini_code[PAR_P1P2];
  double mach = g_inputParam[PAR_MACH] * ini_code[PAR_MACH];

  /* Get downstream mu from T -- not programmed yet
   * need to create new fitting function. Or iterate. */
  static int once1 = 0;
  static double mu2 = 1.3;
    if (!once1) {
      double mu2_old = 0.6;
      v[VX1] = v[VX2] = v[VX3] = 0.0;
      v[TRC] = 1.0;
      v[RHO] = 1.;
      while (fabs(mu2 - mu2_old) > 0.01) {
          v[PRS] = te2 / mu2;
          mu2_old = mu2;
          mu2 = MeanMolecularWeight(v);
      }
      once1 = 1;
  }

  /* Calculate pressures */
  double prs1, prs2, vsnd;
  prs2 = rho2 * te2  / mu2;
  prs1 = p1p2 * prs2 / (g_gamma * mach * mach + 1.);
  vsnd = sqrt(g_gamma * prs1 / rho1);

  /* Calculate p* and u* so that we can transform by -u* and be in the frame of the contact discontinuity */

  static int once2 = 0;
  static double prs_star, u_star;
    if (!once2) {
      double q_l[] = {rho1, vsnd * mach, prs1};
      double q_r[] = {rho2, 0, prs2};

      prs_star = calc_p_star(q_l, q_r, g_gamma);
      u_star = calc_u_star(prs_star, q_l, q_r, g_gamma);

      once2 = 1;
  }

  /* Fill primitives array */

  /* Cloud region */
  if (x1 > 0.) {
    v[RHO] = rho2;
    v[VX1] = -u_star;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
#if HAVE_ENERGY
    v[PRS] = prs2;
#endif
    v[TRC] = 1.0;
  }

  /* Hot, diffuse plasma region */
  else {
    v[RHO] = rho1;
    v[VX1] = vsnd * mach - u_star;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
#if HAVE_ENERGY
    v[PRS] = prs1;
#endif
    v[TRC] = 0.0;
  }

  g_minCoolingTemp = 100.;

  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

    /* Find the biggest pressure gradient - that should be location of contact discontinuity */

//    double * grad_rho = ARRAY_1D(NX1, double);
//    double * x1 = grid->x[IDIR];
//    double *rho = d->Vc[RHO][KBEG][JBEG];
//    int k, j, i;
//    int imax;
//    double v_shift, v_sound;
//
//    GradFromArray(rho, x1, grad_rho, NX1);
//    imax = ArgMaxArray(grad_rho, NX1);
//
//    /* Perform a Gallilean transform by the amount of the velocity at imax */
//    v_shift = d->Vc[VX1][0][0][imax];
//    v_sound = sqrt(g_gamma * d->Vc[PRS][0][0][imax] / d->Vc[RHO][0][0][imax]);
//    TOT_LOOP(k, j, i) {
//        d->Vc[VX1][k][j][i] -= (v_shift + v_sound);
//    }
//
//    RBox box;
//    RBoxDefine (0, NX1_TOT, 0, NX2_TOT, 0, NX3_TOT, CENTER, &box);
//    PrimToCons3D(d->Vc, d->Uc, &box);

//    printf("Shifting by v_shift + v_sound = %g\n", (v_shift + v_sound) * vn.v_norm / 1.e5);

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
