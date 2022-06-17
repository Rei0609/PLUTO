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
    double radius;
    double area;
    double rho;
    double prs;
} AGN;
AGN agn;

typedef struct {
    double nhot;
    double Thot;
    double nwarm;
    double Twarm;
    double rwarm;
    double sigma;
    double Tcrit;
    double te;
} ISM;
ISM ism;

//  1 time scale = 3.*10^3 year

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
    //　AGN initial condition //
    double mach, power, speed, area;
    double rho0, prs0, r_in, q;
    static int agn_once = 0;

    if (! agn_once) {
        /* Input parameters */
        mach = g_inputParam[PAR_MACH];
        power = g_inputParam[PAR_POWER];
        speed = g_inputParam[PAR_SPEED];

        /* Normalize to code units */
        power /= UNIT_VELOCITY * UNIT_VELOCITY * UNIT_VELOCITY *
                 UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH;

        r_in = g_inputParam[PAR_RADIUS];

        /* Flux through surface */
        area = CONST_PI * r_in * r_in;

        /* Primitive variables from input parameters */
        q = (g_gamma - 1.) * mach * mach;
        q = 1. + 2. / q;

        /* Init density and pressure of AGN */
        rho0 = 2. * power / (speed * speed * speed * area * q);
        prs0 = (2. * power / speed - rho0 * speed * speed * area) *
               (g_gamma - 1.) / (2. * g_gamma * area);

        agn.power = power;
        agn.mach = mach;
        agn.speed = speed;
        agn.radius = r_in;
        agn.area = area;
        agn.rho = rho0;
        agn.prs = prs0;

        agn_once = 1;
    }

    // calculation area initial condition //
    static int ism_once = 0;

    if (! ism_once) {

        ism.rwarm = g_inputParam[PAR_RADIUS_W];
        ism.nwarm = g_inputParam[PAR_RHOW];
        ism.nhot = g_inputParam[PAR_RHOH];
        ism.Thot = g_inputParam[PAR_THOT] / KELVIN;
        ism.sigma = g_inputParam[PAR_SIGMA];
        ism.Tcrit = g_inputParam[PAR_TCRIT] / KELVIN;

        ism_once = 1;
    }


  v[RHO] = ism.nhot;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = v[RHO] * ism.Thot /0.6063;
  v[TRC] = 0.0;
  v[TRC+1] = 0.0;

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
    int i, j, k;
    int id, size[3];
    double r, r2, r3, f_exp, f1, f2, f, nc, trw, T_w;
    double r_f, focus, dis, r_1, r_2;
    long int offset, offset1;
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
//    double *x1r = grid->xr[IDIR];
//    double *x2r = grid->xr[JDIR];
//    double *x3r = grid->xr[KDIR];

    /* Interpolate density */
    id = InputDataOpen("input-rho_128_12.flt", "grid_in_128.out", " ", 0, CENTER);


    /* Smooth cloud region */
    DOM_LOOP(k, j, i) {

                /* Get fractal cube */
                nc = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
                /* Apodize with tapered profile */
                r = DIM_EXPAND(x1[i] * x1[i], + x2[j] * x2[j], + x3[k] * x3[k]);
                r = sqrt(r);
//              f_exp = exp((r - ism.rwarm) / ism.sigma);
//              f1 = 1. - f_exp;
//              f2 = f_exp;
//              nc *= f1 * ism.nhot + f2 * ism.nwarm;
                f = 1 + exp((r - ism.rwarm) / ism.sigma);
                f = (log10(ism.nwarm) - log10(ism.nhot)) * 1/f + log10(ism.nhot);
                nc *= pow(10, f);
//                f = (ism.nwarm - ism.nhot) * 1/f + ism.nhot;
//                nc *= f;

                /* Apply critical temperature criterion */
                T_w = ism.nhot / nc * ism.Thot;
                nc = T_w > ism.Tcrit ? ism.nhot : nc;
                trw = T_w > ism.Tcrit ? 0 : 1;

        #if START_MODE == START_MODE_CLOUD

                d->Vc[RHO][k][j][i] = nc;
                d->Vc[PRS][k][j][i] = nc * ism.Thot / 0.6063;
                d->Vc[TRC][k][j][i] = 0;
                d->Vc[TRC + 1][k][j][i] = 1;

        #elif START_MODE == START_MODE_AGN
             if (agn.radius <= 0.1) {

                r3 = sqrt(x1[i] * x1[i] + x2[j] * x2[j] + x3[k] * x3[k]);

                d->Vc[RHO][k][j][i] = nc;
                d->Vc[PRS][k][j][i] = ism.nhot * ism.Thot / 0.6063;
                d->Vc[TRC][k][j][i] = 0;
                d->Vc[TRC + 1][k][j][i] = trw;

                if (r3 <= agn.radius) {
                    d->Vc[RHO][k][j][i] = agn.rho;
                    d->Vc[PRS][k][j][i] = agn.prs;
                    d->Vc[VX3][k][j][i] = agn.speed;
                    d->Vc[TRC][k][j][i] = 1;
                    d->Vc[TRC + 1][k][j][i] = 0;
                }
            } else {

                r2 = sqrt(x1[i] * x1[i] + x2[j] * x2[j]);

                /* r_f is focus spot of ellipsoid */
                r_f = 0.1;

                /* Ellipsoidal shell generation */
                focus = sqrt(agn.radius * agn.radius - r_f * r_f);
                r_1 = r2 + focus;
                r_1 = sqrt(x3[k] * x3[k] + r_1 * r_1);
                r_2 = r2 - focus;
                r_2 = sqrt(x3[k] * x3[k] + r_2 * r_2);
                dis = r_1 + r_2;

                d->Vc[RHO][k][j][i] = ism.nhot;
                d->Vc[PRS][k][j][i] = ism.nhot * ism.Thot / 0.6063;
                d->Vc[TRC][k][j][i] = 0;
                d->Vc[TRC + 1][k][j][i] = 0;

                if (dis <= 2*agn.radius) {
                    d->Vc[RHO][k][j][i] = agn.rho;
                    d->Vc[PRS][k][j][i] = agn.prs;
                    d->Vc[VX3][k][j][i] = agn.speed;
                    d->Vc[TRC][k][j][i] = 1;
                    d->Vc[TRC + 1][k][j][i] = 0;
                }

            }



        #elif START_MODE == START_MODE_HALO

            if (agn.radius <= 0.1) {

                r3 = sqrt(x1[i] * x1[i] + x2[j] * x2[j] + x3[k] * x3[k]);

                d->Vc[RHO][k][j][i] = ism.nhot;
                d->Vc[PRS][k][j][i] = ism.nhot * ism.Thot / 0.6063;
                d->Vc[TRC][k][j][i] = 0;
                d->Vc[TRC + 1][k][j][i] = 0;

                if (r3 <= agn.radius) {
                    d->Vc[RHO][k][j][i] = agn.rho;
                    d->Vc[PRS][k][j][i] = agn.prs;
                    d->Vc[VX3][k][j][i] = agn.speed;
                    d->Vc[TRC][k][j][i] = 1;
                    d->Vc[TRC + 1][k][j][i] = 0;
                }
            } else {

                r2 = sqrt(x1[i] * x1[i] + x2[j] * x2[j]);

                /* r_f is focus spot of ellipsoid */
                r_f = 0.1;

                /* Ellipsoidal shell generation */
                focus = sqrt(agn.radius * agn.radius - r_f * r_f);
                r_1 = r2 + focus;
                r_1 = sqrt(x3[k] * x3[k] + r_1 * r_1);
                r_2 = r2 - focus;
                r_2 = sqrt(x3[k] * x3[k] + r_2 * r_2);
                dis = r_1 + r_2;

                d->Vc[RHO][k][j][i] = ism.nhot;
                d->Vc[PRS][k][j][i] = ism.nhot * ism.Thot / 0.6063;
                d->Vc[TRC][k][j][i] = 0;
                d->Vc[TRC + 1][k][j][i] = 0;

                if (dis <= 2*agn.radius) {
                    d->Vc[RHO][k][j][i] = agn.rho;
                    d->Vc[PRS][k][j][i] = agn.prs;
                    d->Vc[VX3][k][j][i] = agn.speed;
                    d->Vc[TRC][k][j][i] = 1;
                    d->Vc[TRC + 1][k][j][i] = 0;
                }

            }



        #endif
    }

    InputDataClose(id);





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
  int   i, j, k, n;
  double  r2, *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];


  if (side == 0) {    /* -- check solution inside domain -- */
      TOT_LOOP(k,j,i){
    }
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
            BOX_LOOP(box,k,j,i) {
#if START_MODE == START_MODE_AGN
                        r2 = sqrt(x1[i] * x1[i] + x2[j] * x2[j]);
                        if (r2 <= agn.radius) {
                            d->Vc[RHO][k][j][i] = agn.rho;
                            d->Vc[PRS][k][j][i] = agn.prs;
                            d->Vc[VX1][k][j][i] = 0;
                            d->Vc[VX2][k][j][i] = 0;
                            d->Vc[VX3][k][j][i] = agn.speed;
                            d->Vc[TRC][k][j][i] = 1.;
                            d->Vc[TRC+1][k][j][i] = 0;
                        } else {
                            NVAR_LOOP(n) d->Vc[n][k][j][i] = d->Vc[n][2 * KBEG - k - 1][j][i];
                            d->Vc[VX3][k][j][i] *= -1;
                        }
#elif START_MODE == START_MODE_CLOUD
                            NVAR_LOOP(n) d->Vc[n][k][j][i] = d->Vc[n][2 * KBEG - k - 1][j][i];
                            d->Vc[VX3][k][j][i] *= -1;
#elif START_MODE == START_MODE_HALO
                        r2 = sqrt(x1[i] * x1[i] + x2[j] * x2[j]);
                        if (r2 <= agn.radius) {
                            d->Vc[RHO][k][j][i] = agn.rho;
                            d->Vc[PRS][k][j][i] = agn.prs;
                            d->Vc[VX1][k][j][i] = 0;
                            d->Vc[VX2][k][j][i] = 0;
                            d->Vc[VX3][k][j][i] = agn.speed;
                            d->Vc[TRC][k][j][i] = 1.;
                            d->Vc[TRC+1][k][j][i] = 0;
                        } else {
                            NVAR_LOOP(n) d->Vc[n][k][j][i] = d->Vc[n][2 * KBEG - k - 1][j][i];
                            d->Vc[VX3][k][j][i] *= -1;
                        }

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
