//
// Created by Alexander Y. Wagner on 10/30/16.
//

#include "shock.h"
#include "pluto.h"

double mach_from_te1(const double spd, const double te, const double mu) {

    return sqrt( -(4 * (g_gamma - 1) * spd * spd ) /
                 ((1 + (g_gamma - 6) * g_gamma) * spd * spd +
                  (1 + g_gamma) * spd * sqrt(-8 * (g_gamma - 1) * te * g_gamma / mu +
                                                     pow(1 + g_gamma, 2) * spd * spd)));
}

double te_jump(const double mach) {

    return (2 * g_gamma * mach * mach - (g_gamma - 1)) * ((g_gamma - 1) * mach * mach + 2) / \
               pow((g_gamma + 1) * mach, 2);

}

double rho_jump(const double mach) {

    return (g_gamma + 1) * mach * mach / ((g_gamma - 1) * mach * mach + 2);

}

double rho_jump_cdl(const double mach, const double t0, const double td) {

    double tr = t0 / td;

    double p = (g_gamma * mach * mach + 1.) * tr;
    double q = g_gamma * mach * mach * tr;

    return p / 2. + sqrt(p * p / 4. - q);
    // return g_gamma * mach * mach; // This is only in the isothermal and strong shock limit

}

double mach(const double u, const double te, const double mu) {
    return u / sqrt(g_gamma * te / mu);
}

ShockState sh;


double f_shock(double x, double rho, double prs, double gamma) {

    double a, b;

    a = 2. / ((gamma + 1.) * rho);
    b = (gamma - 1.) / (gamma + 1.) * prs;

    return (x - prs) * sqrt(a / (x + b));

}

double f_rarefaction(double x, double rho, double prs, double gamma) {

    double a;

    a = sqrt(gamma * prs / rho);

    return (2. * a) / (gamma - 1.) * (pow(x / prs, (gamma - 1.) / (2. * gamma)) - 1.);
}


double p_equation(double x, double * q_l, double * q_r, double gamma) {

    double rho_l = q_l[0];
    double u_l = q_l[1];
    double prs_l =  q_l[2];
    double rho_r = q_r[0];
    double u_r = q_r[1];
    double prs_r = q_r[2];

    double f_l, f_r;

    if (x > prs_l) f_l = f_shock(x, rho_l, prs_l, gamma);
    else f_rarefaction(x, rho_l, prs_l, gamma);

    if (x > prs_r) f_r = f_shock(x, rho_r, prs_r, gamma);
    else f_rarefaction(x, rho_r, prs_r, gamma);

    return f_l + f_r + u_r - u_l;
}


double secant(double eps, double * q_l, double * q_r, double gamma) {

    double x1, x2, xm, x0;

//    x1 = q_l[2] + q_l[0] * q_l[1] * q_l[1];
//    x2 = q_r[2];
    double xx;
    xx = q_l[2] + q_l[0] * q_l[1] * q_l[1];
    x1 = 3. * xx;
    x2 = xx;

    double p_eq_1, p_eq_2;
    p_eq_1 = p_equation(x1, q_l, q_r, gamma);
    p_eq_2 = p_equation(x2, q_l, q_r, gamma);

    if (p_eq_1 * p_eq_2 < 0) {

        int n = 0;
        do {
            /* Calculate the intermediate value */
            p_eq_1 = p_equation(x1, q_l, q_r, gamma);
            p_eq_2 = p_equation(x2, q_l, q_r, gamma);
            x0 = (x1 * p_eq_2 - x2 * p_eq_1) / (p_eq_2 - p_eq_1);

            /* Update interval and number of iterartions */
            x1 = x2;
            x2 = x0;
            n++;

            /* If x0 is the root of equation then break the loop */

            p_eq_1 = p_equation(x1, q_l, q_r, gamma);
            p_eq_2 = p_equation(x2, q_l, q_r, gamma);

            if (p_eq_1 * p_eq_2 == 0) break;

            xm = (x1 * p_eq_2 - x2 * p_eq_1) / (p_eq_2 - p_eq_1);

        } while (fabs(xm - x0) > eps);

    }
    else {

        printf("Cannot find a root in the given interval");
        QUIT_PLUTO(1)

    }

    return xm;
}


double calc_p_star(double * q_l, double * q_r, double gamma) {


//    double guess = (prs_r + prs_l) / 2.;

    double prs_star = secant(1.e-15, q_l, q_r, gamma);

    return prs_star;

}


/* ********************************************************************* */
double calc_u_star(double prs_star, double * q_l, double * q_r, double gamma) {
/*!
 * Calculates the velocity in the intermediate region
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [in] prs_star  pressure in star region
 * \param [in] q_l       distant left state vector
 * \param [in] q_r       distant right state vector
 * \param [in] gamma     adiabatic index
*********************************************************************** */

    double rho_l = q_l[0];
    double u_l = q_l[1];
    double prs_l =  q_l[2];
    double rho_r = q_r[0];
    double u_r = q_r[1];
    double prs_r = q_r[2];

    int wave_type_l, wave_type_r;
    double f_l, f_r;

    if (prs_star > prs_l) {
        wave_type_l = 0;
        f_l = f_shock(prs_star, rho_l, prs_l, gamma);
    }
    else {
        wave_type_l = 1;
        f_l = f_rarefaction(prs_star, rho_l, prs_l, gamma);
    }

    if (prs_star > prs_r) {
        wave_type_r = 0;
        f_r = f_shock(prs_star, rho_r, prs_r, gamma);
    }
    else {
        wave_type_r = 1;
        f_r = f_rarefaction(prs_star, rho_r, prs_r, gamma);
    }

    return (u_l + u_r) / 2. + (f_r - f_l) / 2.;
//    return wave_type_l;
//    return wave_type_r;
}
