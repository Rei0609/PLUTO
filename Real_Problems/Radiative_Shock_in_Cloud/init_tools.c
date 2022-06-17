/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains functions to assist problem initialization in init.c.

  The init_tools.c file contains helper functions used mainly in init.c. This is to keep init.c relatively clean, with only functions that were originally defined. 

  \author AYW 
  \date   2012-12-07 14:41 JST
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"
#include "init_tools.h"
#include "shock.h"

/* Global struct and arrays for normalization */
VarNorm vn;
double ini_cgs[USER_DEF_PARAMETERS];
double ini_code[USER_DEF_PARAMETERS];


/* Functions */

/* ************************************************************** */
void PrintInitData01(){
/*
 * Print some additional data during initialization
 *
 **************************************************************** */

  double rho1 = g_inputParam[PAR_RHO1] * ini_code[PAR_RHO1];
  double rho2 = g_inputParam[PAR_RHO2] * ini_code[PAR_RHO2];
  double te2 = g_inputParam[PAR_TE2] * ini_code[PAR_TE2];
  double p1p2 = g_inputParam[PAR_P1P2] * ini_code[PAR_P1P2];


    print("\n");
    print("> ShockTube:\n");
    print("\n");
//    print("region      %14s  %14s  %14s\n", "rho", "u", "te");
//    print("upstream    %14g  %14g  %14g\n", g_inputParam[PAR_RHO1], 0., te1);
//    print("cloud       %14g  %14g  %14g\n", g_inputParam[PAR_RHO2], 0., te2);
//    print("\n");
}


/* ************************************************************** */
void SetBaseNormalization() {
/*
 * Sets initializes VarNorm struct with derived normalizations
 * Gives cgs units upon multiplication form code units.
 *
 * Note that pot_norm is also the (mass) specific energy density [erg / g],
 * eint_norm is the energy [erg], and pres_norm is energy density [erg / cm^-3]
 *
 **************************************************************** */

    vn.l_norm = UNIT_LENGTH;
    vn.dens_norm = UNIT_DENSITY;
    vn.v_norm = UNIT_VELOCITY;
    vn.temp_norm = KELVIN;

    /* Derived normalizations */
    vn.t_norm = vn.l_norm / vn.v_norm;
    vn.area_norm = vn.l_norm * vn.l_norm;
    vn.pres_norm = vn.dens_norm * vn.v_norm * vn.v_norm;
    vn.power_norm = vn.pres_norm * vn.v_norm * vn.area_norm;
    vn.eflux_norm = vn.pres_norm * vn.v_norm;
    vn.eint_norm = vn.pres_norm * vn.l_norm * vn.l_norm * vn.l_norm;
    vn.mdot_norm = vn.dens_norm * pow(vn.l_norm, 3) / vn.t_norm;
    vn.newton_norm = 1. / (vn.t_norm * vn.t_norm * vn.dens_norm);
    vn.pot_norm = vn.v_norm * vn.v_norm;
    vn.acc_norm = vn.v_norm / vn.t_norm;
    vn.n_norm = 1. / (vn.l_norm * vn.l_norm * vn.l_norm);
    vn.m_norm = vn.dens_norm * vn.l_norm * vn.l_norm * vn.l_norm;

    print("> Base normalization initialized.\n\n");

}


/* ************************************************************** */
void SetIniNormalization() {
/*
 * Sets noramlizations for ini input variables.
 * ini_cgs converts ini parameters to cgs upon multiplication,
 * whereas ini_code converts ini parameters to code units.
 *
 **************************************************************** */

    double year, degrad;
    year = CONST_ly / CONST_c;
    degrad = CONST_PI / 180.;

    ini_cgs[PAR_RHO1] = vn.dens_norm;
    ini_cgs[PAR_RHO2] = vn.dens_norm;
    ini_cgs[PAR_TE2] = vn.temp_norm;
    ini_cgs[PAR_P1P2] = 1.;
    ini_cgs[PAR_MACH] = 1.;

    ini_code[PAR_RHO1] = ini_cgs[PAR_RHO1] / vn.dens_norm;
    ini_code[PAR_RHO2] = ini_cgs[PAR_RHO2] / vn.dens_norm;
    ini_code[PAR_TE2] = ini_cgs[PAR_TE2] / vn.temp_norm;
    ini_code[PAR_P1P2] = 1.;
    ini_code[PAR_MACH] = 1.;

    print("> Ini parameter normalization array initialized.\n\n");

    return;

}



/* ************************************************ */
void PrintBaseNormalizations() {
    /*!
     * The function prints out grid structure members and
     * is useful for parallel debugging.
     *
     ************************************************** */

    print("vn.l_norm      = %16e \n", vn.l_norm );
    print("vn.dens_norm   = %16e \n", vn.dens_norm );
    print("vn.v_norm      = %16e \n", vn.v_norm );
    print("vn.temp_norm   = %16e \n", vn.temp_norm );
    print("vn.t_norm      = %16e \n", vn.t_norm );
    print("vn.area_norm   = %16e \n", vn.area_norm );
    print("vn.pres_norm   = %16e \n", vn.pres_norm );
    print("vn.power_norm  = %16e \n", vn.power_norm );
    print("vn.eflux_norm  = %16e \n", vn.eflux_norm );
    print("vn.eint_norm   = %16e \n", vn.eint_norm );
    print("vn.mdot_norm   = %16e \n", vn.mdot_norm );
    print("vn.newton_norm = %16e \n", vn.newton_norm );
    print("vn.pot_norm    = %16e \n", vn.pot_norm );
    print("vn.acc_norm    = %16e \n", vn.acc_norm );
    print("vn.n_norm      = %16e \n", vn.n_norm );
    print("vn.m_norm      = %16e \n", vn.m_norm );
    print("\n");
    print("\n");
    print("\n");

}


/* ************************************************ */
void DxFromXArray(double *x, double *dx, int nx) {
/*!
 * x     input: array of positions
 * dx    output: array of cell widths
 * nx    number of elements in x
 *
 * This function returns an array of cell widths from
 * an array of monotonically increasing 1-D positions.
 * The values of dx at 0 and nx - 1 are properly treated.
 *
 ************************************************** */
    int i;
    for (i = 0; i < nx; i++) {
        double xl = i == 0 ? x[0] : x[i - 1];
        double xr = i == nx - 1 ? x[nx - 1] : x[i + 1];
        dx[i] = (i == 0 || i == nx - 1) ? xr - xl : 0.5 * (xr - xl);
    }
}


/* ************************************************ */
void GradFromArray(double *v, double *x, double *grad, int nx) {
/*!
 * v     input: array of values
 * x     input: array of coordinates where values v exist
 * grad  output: array of cell widths
 * nx    number of elements in v
 *
 * This function returns an array of gradients that are the
 * average of a left and right gradient with weight going to
 * *smaller* gradients. From a function AYW learnt from
 * Sam A. E. G. Falle.
 *
 *     grad = (r r l + r l l) / (r * r + l * l);
 *
 * where l and r are the left and right gradients
 *
 ************************************************** */
    int i;
    for (i = 0; i < nx; i++) {
        double gl = i == 0 ? 0. : (v[i] - v[i - 1]) / (x[i] - x[i - 1]);
        double gr = i == nx - 1 ? 0. : (v[i + 1] - v[i]) / (x[i + 1] - x[i]);
        if (fabs(gr * gr + gl * gl) > 1.e-27) {
            grad[i] = (gr * gr * gl + gr * gl * gl) / (gr * gr + gl * gl);
        }
        else {
            grad[i] = (gr + gl) / 2.;
        }
    }
}


/* ************************************************ */
int ArgMaxArray(const double * v, int nx) {
/*!
 * v     input: array of values
 * x     input: array of coordinates where values v exist
 *
 * Return index of maximum value from a 1D array
 *
 ************************************************** */

    double v_old = -1.e30;
    int i;
    int imax = 0;
    for (i = 0; i < nx; i++) {
        if (v[i] > v_old) {
            v_old = v[i];
            imax = i;
        }
    }

    return i;

}


void TransposeArray(double** src, double** dst, int n, int m) {
    int i, j;
    for(i = 0; i < n; ++i)
        for(j = 0; j < m; ++j)
            dst[j][i] = src[i][j];
}

