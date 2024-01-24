
/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>

#include "func_spiral_2b.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int N_PARAMS = 23;
const char PARAM_LABELS[][30] = {"PA", "ell", "r_0", "phi_0", "r_break_1", "phi_break_1", "r_break_2", "phi_break_2", "r_end", "phi_end", "mu_a_2", "mu_b_2", "mu_c_2",
                                 "I_0", "part_growth", "h_s", "part_cutoff",
                                 "width", "w_asymm", "n_out", "n_in", "gamma_out", "gamma_in"};

const char FUNCTION_NAME[] = "Spiral arm function with 2 breaks";
const double pi = 3.14159265358979323846;
const double pi2 = 6.28318530718;
const double DEG2RAD = 0.017453292519943295;

const char SpiralArm2b::className[] = "SpiralArm2b";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

SpiralArm2b::SpiralArm2b() {
    string paramName;

    nParams = N_PARAMS;
    functionName = FUNCTION_NAME;
    shortFunctionName = className;

    // Set up the vector of parameter labels
    for (int i = 0; i < nParams; i++) {
        paramName = PARAM_LABELS[i];
        parameterLabels.push_back(paramName);
    }

    doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void SpiralArm2b::Setup(double params[], int offsetIndex, double xc, double yc) {
    x0 = xc;
    y0 = yc;
    PA = (params[0 + offsetIndex] + 90.0) * DEG2RAD;
    ell = params[1 + offsetIndex];
    r_0 = params[2 + offsetIndex];
    phi_0 = params[3 + offsetIndex] * DEG2RAD;
    r_break_1 = params[4 + offsetIndex];
    phi_break_1 = params[5 + offsetIndex] * DEG2RAD;
    r_break_2 = params[6 + offsetIndex];
    phi_break_2 = params[7 + offsetIndex] * DEG2RAD;
    r_end = params[8 + offsetIndex];
    phi_end = params[9 + offsetIndex] * DEG2RAD;
    mu_a_2 = params[10 + offsetIndex];
    mu_b_2 = params[11 + offsetIndex];
    mu_c_2 = params[12 + offsetIndex];
    I_0 = params[13 + offsetIndex];
    part_growth = params[14 + offsetIndex];
    h_s = params[15 + offsetIndex];
    part_cutoff = params[16 + offsetIndex];
    width = params[17 + offsetIndex];
    w_asymm = params[18 + offsetIndex];
    n_out = params[19 + offsetIndex];
    n_in = params[20 + offsetIndex];
    gamma_out = params[21 + offsetIndex];
    gamma_in = params[22 + offsetIndex];
    
    // pre-compute useful things for this round of invoking the function
    q = 1.0 - ell;
    cosPA = cos(PA);
    sinPA = sin(PA);
    
    is_clockwise = phi_0 - phi_end;
    psi_end = abs(is_clockwise);

    if (is_clockwise > 0) {
        psi_break_1 = phi_0 - phi_break_1;
        psi_break_2 = phi_0 - phi_break_2;
        if (phi_0 >= 0) {
            m_phi_0 = pi - phi_0;
        } else {
            m_phi_0 = -phi_0 - pi;
        }
    } else {
        psi_break_1 = phi_break_1 - phi_0;
        psi_break_2 = phi_break_2 - phi_0;
        if (phi_0 <= pi) {
            m_phi_0 = phi_0;
        } else {
            m_phi_0 = phi_0 - pi2;
        }
    }
    
    psi_growth = part_growth * psi_end;
    psi_cutoff = part_cutoff * psi_end;
    
    mu_a_1 = log(r_break_1 / r_0);
    mu_b_1 = log(r_break_2 / r_break_1);
    mu_c_1 = log(r_end / r_break_2);
    
    m_a_1 = mu_a_1 - mu_a_2;
    m_a_2 = mu_a_2;
    m_b_1 = mu_b_1 - mu_b_2;
    m_b_2 = mu_b_2;
    m_c_1 = mu_c_1 - mu_c_2;
    m_c_2 = mu_c_2;
    
    w_out = width * (1 + w_asymm) / 2;
    w_in = width * (1 - w_asymm) / 2;
    
    n_out_inv = 1.0 / n_out;
    n_in_inv = 1.0 / n_in;
    bn = log(2);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double SpiralArm2b::GetValue(double x, double y) {
    double xc = x - x0;
    double yc = y - y0;

    double xp = xc * cosPA + yc * sinPA;
    double yp = -xc * sinPA + yc * cosPA;

    double x_d0 = xp;
    double y_d0 = yp / q;

    double r = sqrt(x_d0 * x_d0 + y_d0 * y_d0);

    if (is_clockwise > 0)
        x_d0 = -x_d0;

    double psi = atan2l(y_d0, x_d0) - m_phi_0;

    if (psi < 0)
        psi += pi2;
    if (psi >= pi2)
        psi -= pi2;

    double psi_in = GetNearestCoordinates(r, psi);
    double psi_out;
    if (psi_in == -1) {
        psi_out = psi;
    } else
        psi_out = psi_in + pi2;

    while (psi_in > psi_out)
        psi_in -= pi2;

    double I = 0;
    I += GetBrightness(psi_in, r);
    I += GetBrightness(psi_out, r);
    return I;
}

/* ----------------------------- OTHER FUNCTIONS -------------------------------- */

double SpiralArm2b::GetBrightness(double psi, double r) {
    if (psi <= 0 || psi >= psi_end)
        return 0;

    double r_spiral = GetRadius(psi);
    double I;
    I = I_0 * exp(-r_spiral / h_s) * GetNormalBrightness(psi, r - r_spiral);

    if (psi < psi_growth) {
        return I * (3 * pow(psi / psi_growth, 2) - 2 * pow(psi / psi_growth, 3));
    } else if (psi > psi_end - psi_cutoff) {
        return I * ((psi_end - psi) / (psi_cutoff));
    } else
        return I;
}


double SpiralArm2b::GetNormalBrightness(double psi, double h) {
    if (psi == 0)
        return 0;
    double w = h < 0 ? w_in : w_out;
    double n_inv = h < 0 ? n_in_inv : n_out_inv;
    double gamma = h < 0 ? gamma_in : gamma_out;
    double loc_w = w * exp(gamma * (psi / psi_end - 0.5));
    return exp(-bn * pow(abs(h) / loc_w, n_inv));
}

double SpiralArm2b::GetRadius(double psi) {
    double psi_norm;
    if (psi >= psi_break_2) {
        psi_norm = (psi - psi_break_2) / (psi_end - psi_break_2);
    	return r_break_2 * exp(m_c_1 * psi_norm + m_c_2 * pow(psi_norm, 2));
    } else if (psi >= psi_break_1) {
        psi_norm = (psi - psi_break_1) / (psi_break_2 - psi_break_1);
    	return r_break_1 * exp(m_b_1 * psi_norm + m_b_2 * pow(psi_norm, 2));
    } else {
        psi_norm = (psi / psi_break_1);
    	return r_0 * exp(m_a_1 * psi_norm + m_a_2 * pow(psi_norm, 2));
    }
}

double SpiralArm2b::GetNearestCoordinates(double r, double psi) {
    double r_2 = GetRadius(psi);
    double psi_1 = -1;
    double psi_2 = psi;
    int i = 0;
    while (r_2 < r && psi_2 < psi_end && i++ < 15) {
        psi_1 = psi_2;
        psi_2 = psi_2 + pi2;
        r_2 = GetRadius(psi_2);
    }

    return psi_1;
}

/* END OF FILE: func_spiral_2b.cpp ------------------------------------ */
