/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute oiii cooling rate

  \authors A. Mignone (mignone@ph.unito.it)\n
           M. Sormani\n

 \b References

  \date   Apr 09, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* AYW -- 2013-01-08 23:05 JST */
#include "init_tools.h"

extern double gCooling_x1, gCooling_x2, gCooling_x3;

/* ***************************************************************** */
void Radiat_frac (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
    int klo, khi, kmid;
    static int ntab;
    double mu, T, Tmid, scrh, dT, prs;
    static double *L_tab, *T_tab, E_cost;
    static double *L_oiii_tab, *L_ha_tab, *L_x1_tab, *L_x2_tab, *L_x3_tab, *L_x4_tab;
    double nH = UNIT_DENSITY/CONST_amu*H_MASS_FRAC/CONST_AH*v[RHO];
    double ne = nH*(1.0 + 0.5*CONST_AZ*FRAC_Z);

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

    if (T_tab == NULL) {
        FILE *fcool;
        printLog(" > Reading cooling fraction table from disk...\n");

        /*  Read in table */
        /* Losses (L)     L_5007        LHalpha       0.1-0.5keV    0.5-1.0keV    1.0-2.0keV    2.0-10.0keV  */
        /* (erg/cm^3/s)   (erg cm^3/s)  (erg cm^3/s)  (fraction)    (fraction)    (fraction)    (fraction)   */
        fcool = fopen("cooltable_frac.dat", "r");
        if (fcool == NULL) {
            printLog ("! Radiat: cooltable_frac.dat could not be found.\n");
            QUIT_PLUTO(1);
        }
        L_tab = ARRAY_1D(20000, double);
        T_tab = ARRAY_1D(20000, double);
        L_oiii_tab = ARRAY_1D(20000, double);
        L_ha_tab = ARRAY_1D(20000, double);
        L_x1_tab = ARRAY_1D(20000, double);
        L_x2_tab = ARRAY_1D(20000, double);
        L_x3_tab = ARRAY_1D(20000, double);
        L_x4_tab = ARRAY_1D(20000, double);

        ntab = 0;
        while (fscanf(fcool, "%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", T_tab + ntab,
                      L_tab + ntab, L_oiii_tab + ntab, L_ha_tab + ntab,
                      L_x1_tab + ntab, L_x2_tab + ntab, L_x3_tab + ntab, L_x4_tab + ntab) != EOF) {
            ntab++;
        }
        E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
        fclose(fcool);
    }

/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

    /* Intead of T, we find the cooling rate as a function of T = p / rho */
    prs = v[0] * (g_gamma - 1.0);
    if (prs < 0.0) {
        prs = g_smallPressure;
        v[0] = prs / (g_gamma - 1.0);
    }
    mu = MeanMolecularWeight(v);
    T = prs / v[RHO] * KELVIN * mu;

    if (T != T) {
    printf ("! Radiat(): Nan found: rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
        QUIT_PLUTO(1);
    }

/*
  if (T < g_minCoolingTemp) {
    rhs[0] = 0.0;
    return;
  }
*/
/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

    klo = 0;
    khi = ntab - 1;

    if (T > T_tab[khi] || T < T_tab[klo]) {
    rhs[0] = 0.0;
    return;
        QUIT_PLUTO(1);
    }

    while (klo != (khi - 1)) {
        kmid = (klo + khi) / 2;
        Tmid = T_tab[kmid];
        if (T <= Tmid) {
            khi = kmid;
        } else if (T > Tmid) {
            klo = kmid;
        }
    }

/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */

  dT        = T_tab[khi] - T_tab[klo];
  scrh      = L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT;
  rhs[0] = -ne*nH*scrh*E_cost;

/* ----------------------------------------------
    Temperature cutoff
   ---------------------------------------------- */

  rhs[0] *= 1.0 - 1.0/cosh( pow( T/g_minCoolingTemp, 12));

/* ----------------------------------------------
    Other quantities
   ---------------------------------------------- */

    scrh      = L_oiii_tab[klo]*(T_tab[khi] - T)/dT + L_oiii_tab[khi]*(T - T_tab[klo])/dT;
    rhs[1] = -ne*nH*scrh*E_cost;

    scrh      = L_ha_tab[klo]*(T_tab[khi] - T)/dT + L_ha_tab[khi]*(T - T_tab[klo])/dT;
    rhs[2] = -ne*nH*scrh*E_cost;

    scrh      = L_x1_tab[klo]*(T_tab[khi] - T)/dT + L_x1_tab[khi]*(T - T_tab[klo])/dT;
    rhs[3] = rhs[0] * scrh;

    scrh      = L_x2_tab[klo]*(T_tab[khi] - T)/dT + L_x2_tab[khi]*(T - T_tab[klo])/dT;
    rhs[4] = rhs[0] * scrh;

    scrh      = L_x3_tab[klo]*(T_tab[khi] - T)/dT + L_x3_tab[khi]*(T - T_tab[klo])/dT;
    rhs[5] = rhs[0] * scrh;

    scrh      = L_x4_tab[klo]*(T_tab[khi] - T)/dT + L_x4_tab[khi]*(T - T_tab[klo])/dT;
    rhs[6] = rhs[0] * scrh;

}
