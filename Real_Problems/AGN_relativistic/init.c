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

static double Profile(double r, int nv, int xn, double r0);

typedef struct {
    double power;
    double mach;
    double speed;
    double radius;
    double area;
    double rho;
    double prs;
    double gamma;
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
    double gamma;
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
    double mach, power, speed, area, lgamma, lgamma2;
    double rho0, prs0, r_in, q, q_rel, q_non, x, chi, gamma_rel, gamma_non, mach_non;
    static int agn_once = 0;

    if (! agn_once) {
        /* Input parameters */
        mach = g_inputParam[PAR_MACH];
        power = g_inputParam[PAR_POWER];
        lgamma = g_inputParam[PAR_GAMMA];
        gamma_rel = 5./3.;
        gamma_non = 5./3.;

        /* Normalize to code units */
        power /= UNIT_VELOCITY * UNIT_VELOCITY * UNIT_VELOCITY *
                 UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH;

        r_in = g_inputParam[PAR_RADIUS];

        /* Flux through surface */
        area = CONST_PI * r_in * r_in;

        /* Primitive variables from input parameters */
        lgamma2 = lgamma * lgamma;
        speed = sqrt(1 - 1 / lgamma2);
        q = (gamma_non - 1.) * mach * mach;
        q = 1. + 2. / q;
        q_rel = gamma_rel / (gamma_rel - 1.);
        q_non = gamma_non / (gamma_non - 1.);

        /* Init parameter of AGN (RHD) */
        chi = mach * mach * (gamma_rel - 1) / (lgamma2 - 1);
        chi = chi + gamma_rel - 2;
        prs0 = lgamma2 * area * speed * ((1 - 1 / lgamma) * chi + 1) * q_rel;
        prs0 = power / prs0;
        rho0 = prs0 * q_rel * chi;

        /*
        prs0 = chi * (lgamma + 1) / lgamma + 1;
        prs0 = power / (area * lgamma2 * speed * q * chi * prs0);
        rho0 = prs0 * chi * q / lgamma2;
         */

        printf("JetRho:%f\n",rho0);
        printf("JetPrs:%f\n",prs0);
        printf("JetSpd:%f\n",speed);
        printf("Chi:%f\n",chi);
        printf("gamma:%f\n",gamma_rel);
        printf("Power:%f\n",power);

#if PHYSICS == HD
        if (g_inputParam[PAR_NRJET] == 0) {
            /* Init pressure of AGN (density matched) */
            prs0 = power / (speed * area) - rho0 * speed * speed / 2;
            prs0 = prs0 / q_non;
        } else if (g_inputParam[PAR_NRJET] == 1) {
            /* Init density of AGN (pressure matched) */
            rho0 = power / (speed * area) - q_non * prs0;
            rho0 = rho0 * 2 / (speed * speed);
        } else if (g_inputParam[PAR_NRJET] == 2) {
            /* Init pressure and density of AGN (Mach number matched) */
            rho0 = 2. * power / (speed * speed * speed * area * q);
            prs0 = (2. * power / speed - rho0 * speed * speed * area) *
                   (g_gamma - 1.) / (2. * g_gamma * area);
        } else {
            /* Init pressure and density of AGN (Chi matched) */
            prs0 = q_rel * speed * area * (speed * speed * chi / 2 + 1);
            prs0 = power / prs0;
            rho0 = q_rel * prs0 * chi;
        }


        mach_non = speed / sqrt(gamma_non * prs0 / rho0);

        printf("Mach_non:%f\n",mach_non);
        printf("JetRho:%f\n",rho0);
        printf("JetPrs:%f\n",prs0);
        printf("JetSpd:%f\n",speed);

#endif

        agn.power = power;
        agn.mach = mach;
        agn.speed = speed;
        agn.radius = r_in;
        agn.area = area;
        agn.rho = rho0;
        agn.prs = prs0;
        agn.gamma = gamma_rel;

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
        ism.gamma = gamma_non;

        ism_once = 1;
    }


  v[RHO] = ism.nhot;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = v[RHO] * ism.Thot /0.6063;
  v[TRC] = 0.0;
  v[TRC+1] = 0.0;

  g_minCoolingTemp = 100;


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
    id = InputDataOpen("input-rho.flt", "grid_in.out", " ", 0, CENTER);


    /* Smooth cloud region */
    DOM_LOOP(k, j, i) {

        trw = 0;
        nc = ism.nhot;

#if CLOUDS == YES
                /* Get fractal cube */
                nc = InputDataInterpolate(id, x1[i], x2[j], x3[k]);

                /* Apodize with tapered profile */
                r = DIM_EXPAND(x1[i] * x1[i], + x2[j] * x2[j], + x3[k] * x3[k]);
                r = sqrt(r);

                f = 1 + exp((r - ism.rwarm) / ism.sigma);
                f = (log10(ism.nwarm) - log10(ism.nhot)) * 1/f + log10(ism.nhot);
                nc *= pow(10, f);

                /* Apply critical temperature criterion */
                T_w = ism.nhot / nc * ism.Thot;
                nc = T_w > ism.Tcrit ? ism.nhot : nc;
                trw = T_w > ism.Tcrit ? 0 : 1;
#endif

        d->Vc[RHO][k][j][i] = nc;
        d->Vc[PRS][k][j][i] = ism.nhot * ism.Thot / 0.6063;
        d->Vc[TRC][k][j][i] = 0;
        d->Vc[TRC + 1][k][j][i] = trw;


        r2 = sqrt(x1[i] * x1[i] + x2[j] * x2[j]);

        /* r_f is focus spot of ellipsoid */
        r_f = 0.02;

        /* Ellipsoidal shell generation */
        focus = sqrt(agn.radius * agn.radius - r_f * r_f);
        r_1 = r2 + focus;
        r_1 = sqrt(x3[k] * x3[k] + r_1 * r_1);
        r_2 = r2 - focus;
        r_2 = sqrt(x3[k] * x3[k] + r_2 * r_2);
        dis = r_1 + r_2;

        if (dis <= 2. * agn.radius) {
            d->Vc[RHO][k][j][i] = agn.rho;
            d->Vc[PRS][k][j][i] = agn.prs;
            d->Vc[VX3][k][j][i] = agn.speed;
            d->Vc[TRC][k][j][i] = 1;
            d->Vc[TRC + 1][k][j][i] = 0;

        }

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
  int   i, j, k, nv;
  double  r2, *x1, *x2, *x3;
  double vwnd[NVAR], vrfl[NVAR];

  /* Region around the wind inlet where strong rarefactions are artificially shielded. Hardcoded. */
  double vrfl_prof;


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


                        r2 = sqrt(x1[i] * x1[i] + x2[j] * x2[j]);

                        /* Wind primitives and reflective BC primitives */
                        vwnd[RHO] = agn.rho;
                        vwnd[PRS] = agn.prs;
                        vwnd[VX1] = 0;
                        vwnd[VX2] = 0;
                        vwnd[VX3] = agn.speed;
                        vwnd[TRC] = 1.;
                        vwnd[TRC+1] = 0;
                        NVAR_LOOP(nv) vrfl[nv] = d->Vc[nv][2 * KBEG - k - 1][j][i];

                        /* Create a region at radius PAR_RFGUARD that guards against rarefactions */
//                        vrfl_prof = Profile(r2, nv, 4, g_inputParam[PAR_RFGUARD]) - 1.;
//                        vrfl[VX3] *= vrfl_prof;

                        /* Apply edge-smoothed wind profile */
                        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vrfl[nv] + (vwnd[nv] - vrfl[nv]) *
                                Profile(r2, nv, g_inputParam[PAR_WPROF_IDX], agn.radius);

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

/* ********************************************************************* */
double Profile(double r, int nv, int xn, double r0)
/*
 *
 *
 *********************************************************************** */
{

    if (nv == RHO) r0 *= 1.1;

#if GEOMETRY == SPHERICAL
    r0 = 5.0/180.0*CONST_PI;
#endif
    return 1.0 / cosh(pow(r / r0, xn));
}

