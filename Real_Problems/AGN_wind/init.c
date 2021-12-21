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


typedef struct {
    double power;
    double mach;
    double speed;
    double area;
    double rho;
    double prs;
} AGN;


AGN agn;

typedef struct {
    double *rho;
    double *prs;
    double *r;
    double *logr;
    double *pq_sin;
    double *pq_shell;
    int nr;
} HOT_HALO;
HOT_HALO hot;

typedef struct {
    double *psi;
    double *g;
    double *r;
    double *logr;
    int nr;
} DOUBLE_ISOTHERMAL;
DOUBLE_ISOTHERMAL di;


void rhs_poisson(double x, double *y, double *f) {

  /* RHS of poisson equaiton for double isothermal from 
   * Sutherland & Bicknell (2007), Eq. 5. 
   * The second order differential equaiton is split into
   * two first order diff. equations, with y[0] being the potential,
   * and y[1] the differential of the potential (the force). */

  double b, h, kappa, l, rate, lok2;

  b = 2 / x;
  kappa = g_inputParam[PAR_KAPPA];
  l = g_inputParam[PAR_LAMBDA];
  lok2 = l / kappa;
  lok2 = lok2 * lok2;
  h = exp(-y[0]) + lok2 * exp(-kappa * kappa * y[0]);
  h = 9 * h;

  f[0] = y[1];
  f[1] = h - b * y[1];

}

double interpolate_1D(double r, double *rarr, double *varr, int nr) {
    double ratio, ratio2, end, beg;
    double dr, lr, rc;
    int i;

    end = rarr[nr-1];
    beg = rarr[0];

    lr = end - beg;
    dr = lr / (nr - 1);
    lr += dr;
    beg -= 0.5 * dr;
    end += 0.5 * dr;

    rc = r - 0.5 * dr;

    ratio = (rc - beg) / lr;
    i = (int) (ratio * nr);

    ratio2 = (r - rarr[i]) / dr;

    return varr[i] + ratio2 * (varr[i+1] - varr[i]);
}

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
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
  
  double te;
  double r, r_in;

  r = x1;
  r_in = g_domBeg[IDIR];

  te = g_inputParam[PAR_THOT] / KELVIN;

  v[RHO] = g_inputParam[PAR_NHOT];
  v[VX1] = 0;
  v[VX2] = 0;
  v[VX3] = 0;
  v[PRS] = v[RHO] * te / 0.6063; 
  v[TRC] = 0.0;


  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

 /********************************************************************* */
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


  double r, q, r_in;
  double mach, power, speed, area;
  double rho0, vx1, vx2, vx3, prs0;

  /* Input parameters */
  mach = g_inputParam[PAR_MACH];
  power = g_inputParam[PAR_POWER];
  speed = g_inputParam[PAR_SPEED];


  /* Normalize to code units */
  power /= UNIT_VELOCITY * UNIT_VELOCITY * UNIT_VELOCITY * 
	  UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH;
  
  //r = x1;
  r_in = g_domBeg[IDIR];

  /* Flux through surface of a sphere */
  area = 4 * CONST_PI * r_in * r_in;

  /* Primitive variables from input parameters */

  q = (g_gamma - 1.) * mach * mach;
  q = 1. + 2. / q;

  rho0 = 2. * power / (speed * speed * speed * area * q);
  prs0= (2. * power / speed - rho0 * speed * speed * area) * 
	    (g_gamma - 1.) / (2. * g_gamma * area);
  //prs = power / speed * (1. - 1. / q) * (g_gamma - 1.) / (g_gamma * area);

  agn.power = power;
  agn.mach = mach;
  agn.speed = speed;
  agn.area = area;
  agn.rho = rho0;
  agn.prs = prs0;

  
  /* set up garactic programs */

  /* Pointer to RHS function */

  void * rhs;
  rhs = rhs_poisson;

 
  /* Domain info (r and dr) for ODE solver */
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  x1 = grid->x[IDIR];
  //x2 = grid->x[JDIR];
  //x3 = grid->x[KDIR];
  
  /* Define rd */
  double rd, rb;
  double sigma_d, l, kappa; 
  sigma_d = sqrt(g_inputParam[PAR_TVIR] * CONST_kB / CONST_mp) / (UNIT_VELOCITY);
  rb = g_inputParam[PAR_RBAR];
  //rd = 4 * CONST_PI * CONST_G * rho_d;
  //rd = 9 * sigma_d * sigma_d / rd;
  //rd = sqrt(rd);
  kappa = g_inputParam[PAR_KAPPA];
  l = g_inputParam[PAR_LAMBDA];
  rd = l * rb; 


  /* Find minimum dx1 */
  double dx1_min = 1.e30;
  double *dx1, *dx2, *dx3;
  dx1 = grid->dx[IDIR];
  ITOT_LOOP(i) dx1_min = MIN(dx1_min, dx1[i]);

  /* Initial condition and solution to potential.
   * y[0] is the potential, and y[1] the force. */
  double nvar, y[2];
  nvar = 2;
  y[0] = 0.;
  y[1] = 0.;

  double tratio, te;
  double a, lambda;
  double b, c, u;

  hot.rho = ARRAY_1D(NX1_TOT, double);
  hot.prs = ARRAY_1D(NX1_TOT, double);
  hot.r = ARRAY_1D(NX1_TOT, double);
  hot.logr = ARRAY_1D(NX1_TOT, double);
  hot.pq_sin = ARRAY_1D(NX1_TOT, double);
  hot.pq_shell = ARRAY_1D(NX1_TOT, double);
  di.psi = ARRAY_1D(NX1_TOT, double);
  di.g = ARRAY_1D(NX1_TOT, double);
  di.r = ARRAY_1D(NX1_TOT, double);
  di.logr = ARRAY_1D(NX1_TOT, double);

  /* Solve the 1D radial Poisson equation to obtain profiles of the 
   * potential, density, and pressure.
   * We also want the values in the boundaries */
  print ("> Starting setup of gravitational potential... \n");  
  ITOT_LOOP(i) {
    ODE_Solve(y, nvar, dx1_min*1.e-3/rd, x1[i]/rd, 0.2*dx1_min/rd, rhs, ODE_RK4);
    //di.psi[i] = y[0] * sigma_d * sigma_d;
    di.g[i] = y[1] * sigma_d * sigma_d;
    di.r[i] = x1[i];
    di.logr[i] = log10(x1[i]);
    tratio = g_inputParam[PAR_TVIR] / g_inputParam[PAR_THOT];
    te = g_inputParam[PAR_THOT] / KELVIN;

    hot.rho[i] = g_inputParam[PAR_NHOT] * exp(-tratio * y[0]);
    hot.prs[i] = hot.rho[i] * te / 0.6063;
    hot.r[i] = x1[i];
    hot.logr[i] = log10(x1[i]);
  }
  hot.nr = i;
  di.nr = i;
  print ("> Done\n");


 
  double shell;
  double rho_r, prs_r, r_id;
  r_id = x1[i];


  /* Fill grid data arrays with calculated profile */
  DOM_LOOP(k, j, i) {

             //rho_r = interpolate_1D(log10(r_id), hot.logr, hot.rho, hot.nr);
//             d->Vc[RHO][k][j][i] = rho_r;
             d->Vc[RHO][k][j][i] = hot.rho[i];

#if PERT_MODE == PERT_MODE_SHELL
              /* Induce a shell-like perturbation in density */

              b = 3.;
              c = 0.1 * g_inputParam[PAR_RBAR];
              c = log(2) / ( c * c);
              u = 1 * g_inputParam[PAR_RBAR];
              shell = b * exp(-(x1[i] - u) * (x1[i] - u) * c);
              d->Vc[RHO][k][j][i] = hot.rho[i] * (1. + shell);

#elif PERT_MODE == PERT_MODE_WAVE
              /* Induce a wave-like perturbation in velocity */

              a = 1. * sqrt(g_gamma * CONST_kB * g_inputParam[PAR_THOT]  / (CONST_amu * 0.6063)) / UNIT_VELOCITY;
              lambda = 2 * g_inputParam[PAR_RBAR];
              d->Vc[VX1][k][j][i] = a * sin(2 * CONST_PI * log10(x1[i]) / lambda);


#endif

              d->Vc[PRS][k][j][i] = hot.prs[i];
              

          }

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
  double  a, w, t;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){ 
	#if START_MODE == START_MODE_WIND
	  d->Vc[RHO][k][j][i] = agn.rho;
	  d->Vc[VX1][k][j][i] = agn.speed;
	  d->Vc[VX2][k][j][i] = 0.;
	  d->Vc[VX3][k][j][i] = 0.;
	  d->Vc[PRS][k][j][i] = agn.prs;
	  d->Vc[TRC][k][j][i] = 1.;
	#elif START_MODE == START_MODE_HALO
	  d->Vc[RHO][k][j][i] = hot.rho[i];
	  d->Vc[VX1][k][j][i] = 0.;
	  d->Vc[VX2][k][j][i] = 0.;
	  d->Vc[VX3][k][j][i] = 0.;
	  d->Vc[PRS][k][j][i] = hot.prs[i];
	  d->Vc[TRC][k][j][i] = 0.;
        #endif
      }
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
      BOX_LOOP(box,k,j,i){ 
	  d->Vc[RHO][k][j][i] = hot.rho[i];
	  d->Vc[VX1][k][j][i] = 0.;
	  d->Vc[VX2][k][j][i] = 0.;
	  d->Vc[VX3][k][j][i] = 0.;
	  d->Vc[PRS][k][j][i] = hot.prs[i];
	  d->Vc[TRC][k][j][i] = 0.;
     #if PERT_MODE == PERT_MODE_TAIL
      a = 1. * sqrt(g_gamma * CONST_kB * g_inputParam[PAR_THOT]  / (CONST_amu * 0.6063)) / UNIT_VELOCITY;
      t = 500 * 3.15569252e10 * UNIT_VELOCITY / UNIT_LENGTH;
      w = 2 * CONST_PI / t;
//      d->Vc[RHO][k][j][i] = hot.rho[i] * (1 + a * sin(w * g_time));
//      d->Vc[PRS][k][j][i] = hot.prs[i] * (1 + a * sin(w * g_time));
      d->Vc[VX1][k][j][i] = a * sin(w * g_time);
     #endif
      }
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
  double ratio, ratio2, end, beg, g_r;
  double dlogx1, llogx1, logx1, logx1c;
  int i; 

  g_r = interpolate_1D(log10(x1), di.logr, di.g, di.nr);

  g[IDIR] = -1.0 * g_r;
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
  double ratio, end, beg;
  double dlogx1, llogx1, logx1, logx1c;
  int i; 

  end = log10(g_domEnd[IDIR]);
  beg = log10(g_domBeg[IDIR]);

  llogx1 = end - beg;
  dlogx1 = llogx1 / NX1;

  logx1 = log10(x1);
  logx1c = logx1 + 0.5 * dlogx1;

  ratio = (logx1c - beg) / llogx1;
  i = ratio * NX1 + IBEG;

 return di.psi[i];
}
#endif
