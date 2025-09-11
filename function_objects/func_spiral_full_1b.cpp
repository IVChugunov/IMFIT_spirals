
/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>

#include "func_spiral_full_1b.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int N_PARAMS = 27;
const char PARAM_LABELS[][20] = {"PA", "ell", "r_0", "phi_0", "r_break_1", "phi_break_1", "r_end", "phi_end", "mu_a_2", "mu_a_3", "mu_b_2", "mu_b_3",
                                 "I_0", "part_growth", "ih_s", "part_cutoff",
                                 "w_zp", "w_i", "S_0", "S_end", "n",
                                 "A_d_1", "phi_d_1", "w_d_1", "A_d_2", "phi_d_2", "w_d_2"};

const char FUNCTION_NAME[] = "Spiral arm function with 1 break and extended parameters";
const double pi = 3.14159265358979323846;
const double pi2 = 6.28318530718;
const double DEG2RAD = 0.017453292519943295;

const char SpiralArmFull1b::className[] = "SpiralArmFull1b";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

SpiralArmFull1b::SpiralArmFull1b() {
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

void SpiralArmFull1b::Setup(double params[], int offsetIndex, double xc, double yc) {
    x0 = xc;
    y0 = yc;
    PA = (params[0 + offsetIndex] + 90.0) * DEG2RAD;
    ell = params[1 + offsetIndex];
    r_0 = params[2 + offsetIndex];
    phi_0 = params[3 + offsetIndex] * DEG2RAD;
    r_break_1 = params[4 + offsetIndex];
    phi_break_1 = params[5 + offsetIndex] * DEG2RAD;
    r_end = params[6 + offsetIndex];
    phi_end = params[7 + offsetIndex] * DEG2RAD;
    mu_a_2 = params[8 + offsetIndex];
    mu_a_3 = params[9 + offsetIndex];
    mu_b_2 = params[10 + offsetIndex];
    mu_b_3 = params[11 + offsetIndex];
    I_0 = params[12 + offsetIndex];
    part_growth = params[13 + offsetIndex];
    ih_s = params[14 + offsetIndex];
    part_cutoff = params[15 + offsetIndex];
    w_zp = params[16 + offsetIndex];
    w_i = params[17 + offsetIndex];
    S_0 = params[18 + offsetIndex];
    S_end = params[19 + offsetIndex];
    n = params[20 + offsetIndex];
    A_d_1 = params[21 + offsetIndex];
    phi_d_1 = params[22 + offsetIndex] * DEG2RAD;
    w_d_1 = params[23 + offsetIndex] * DEG2RAD;
    A_d_2 = params[24 + offsetIndex];
    phi_d_2 = params[25 + offsetIndex] * DEG2RAD;
    w_d_2 = params[26 + offsetIndex] * DEG2RAD;
    
    // pre-compute useful things for this round of invoking the function
    q = 1.0 - ell;
    cosPA = cos(PA);
    sinPA = sin(PA);
    
    is_clockwise = phi_0 - phi_end;
    psi_end = abs(is_clockwise);

    if (is_clockwise > 0) {
        psi_break_1 = phi_0 - phi_break_1;
        psi_d_1 = phi_0 - phi_d_1;
        psi_d_2 = phi_0 - phi_d_2;
        if (phi_0 >= 0) {
            m_phi_0 = pi - phi_0;
        } else {
            m_phi_0 = -phi_0 - pi;
        }
    } else {
        psi_break_1 = phi_break_1 - phi_0;
        psi_d_1 = phi_d_1 - phi_0;
        psi_d_2 = phi_d_2 - phi_0;
        if (phi_0 <= pi) {
            m_phi_0 = phi_0;
        } else {
            m_phi_0 = phi_0 - pi2;
        }
    }
    
    psi_growth = part_growth * psi_end;
    psi_cutoff = part_cutoff * psi_end;
    
    mu_a_1 = log(r_break_1 / r_0);
    mu_b_1 = log(r_end / r_break_1);
    
    m_a_1 = mu_a_1 - mu_a_2 + mu_a_3;
    m_a_2 = mu_a_2 - 3 * mu_a_3;
    m_a_3 = 2 * mu_a_3;
    m_b_1 = mu_b_1 - mu_b_2 + mu_b_3;
    m_b_2 = mu_b_2 - 3 * mu_b_3;
    m_b_3 = 2 * mu_b_3;
    
    S_i = (S_end - S_0) / (r_end - r_0);
    S_zp = S_0 - (r_0 * S_i);

    n_inv = 1.0 / n;
    bn = log(2);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double SpiralArmFull1b::GetValue(double x, double y) {
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
    I += GetBrightness(r, psi_in);
    I += GetBrightness(r, psi_out);
    return I;
}

/* ----------------------------- OTHER FUNCTIONS -------------------------------- */

double SpiralArmFull1b::GetBrightness(double r, double psi) {
    if (psi <= 0 || psi >= psi_end)
        return 0;

    double r_spiral = GetRadius(psi);
    double I;
    I = I_0 * GetParallelBrightness(r_spiral, psi) * GetNormalBrightness(r_spiral, r - r_spiral, psi);

    return I;
}

double SpiralArmFull1b::GetParallelBrightness(double r_spiral, double psi) {
    double I_par = exp(-r_spiral * ih_s);

    auto smooth = [](double t) {
        t = std::max(0.0, std::min(1.0, t));
        return (3.0 * pow(t, 2) - 2.0 * pow(t, 3));
    };

    if (psi < psi_growth) {
        double t = psi / psi_growth;
        I_par = I_par * smooth(t);
    }
    if (psi > psi_end - psi_cutoff) {
        double t = (psi_end - psi) / psi_cutoff;
        I_par = I_par * smooth(t);
    }

    if (A_d_1 != 0)
        I_par = I_par / (1 + A_d_1 * exp(-pow(psi - psi_d_1,  2) / (2 * pow(w_d_1, 2))));
    if (A_d_2 != 0)
        I_par = I_par / (1 + A_d_2 * exp(-pow(psi - psi_d_2,  2) / (2 * pow(w_d_2, 2))));

    return I_par;
}

double SpiralArmFull1b::GetNormalBrightness(double r_spiral, double rho, double psi) {
    double S_part = rho < 0 ? (1 - (S_i * r_spiral + S_zp)) / 2 : (1 + (S_i * r_spiral + S_zp)) / 2;
    double loc_w = abs(((w_i * r_spiral) + w_zp) * S_part);
    return exp(-bn * pow(abs(rho) / loc_w, n_inv));
}

double SpiralArmFull1b::GetNearestCoordinates(double r, double psi) {
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

double SpiralArmFull1b::GetRadius(double psi) {
    double psi_norm;
    if (psi >= psi_break_1) {
        psi_norm = (psi - psi_break_1) / (psi_end - psi_break_1);
    	return r_break_1 * exp(m_b_1 * psi_norm + m_b_2 * pow(psi_norm, 2) + m_b_3 * pow(psi_norm, 3));
    } else {
        psi_norm = (psi / psi_break_1);
    	return r_0 * exp(m_a_1 * psi_norm + m_a_2 * pow(psi_norm, 2) + m_a_3 * pow(psi_norm, 3));
    }
}

/* END OF FILE: func_spiral_full_1b.cpp ------------------------------------ */
