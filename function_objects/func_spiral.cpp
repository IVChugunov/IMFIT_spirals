
/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>

#include "func_spiral.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int N_PARAMS = 15;
const char PARAM_LABELS[][20] = {"PA", "ell", "r_0", "phi_0", "r_end", "phi_end", "mu_a_2", "mu_a_3", "mu_a_4",
                                 "I_0", "part_growth", "ih_s", "part_cutoff",
                                 "w_zp", "w_i"};

const char FUNCTION_NAME[] = "Normal spiral arm function with 0 breaks";
const double pi = 3.14159265358979323846;
const double pi2 = 6.28318530718;
const double DEG2RAD = 0.017453292519943295;

const char SpiralArm::className[] = "SpiralArm";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

SpiralArm::SpiralArm() {
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

void SpiralArm::Setup(double params[], int offsetIndex, double xc, double yc) {
    x0 = xc;
    y0 = yc;
    PA = (params[0 + offsetIndex] + 90.0) * DEG2RAD;
    ell = params[1 + offsetIndex];
    r_0 = params[2 + offsetIndex];
    phi_0 = params[3 + offsetIndex] * DEG2RAD;
    r_end = params[4 + offsetIndex];
    phi_end = params[5 + offsetIndex] * DEG2RAD;
    mu_a_2 = params[6 + offsetIndex];
    mu_a_3 = params[7 + offsetIndex];
    mu_a_4 = params[8 + offsetIndex];
    I_0 = params[9 + offsetIndex];
    part_growth = params[10 + offsetIndex];
    ih_s = params[11 + offsetIndex];
    part_cutoff = params[12 + offsetIndex];
    w_zp = params[13 + offsetIndex];
    w_i = params[14 + offsetIndex];
    
    // pre-compute useful things for this round of invoking the function
    q = 1.0 - ell;
    cosPA = cos(PA);
    sinPA = sin(PA);
    
    is_clockwise = phi_0 - phi_end;
    psi_end = abs(is_clockwise);

    if (is_clockwise > 0) {
        m_phi_0 = pi - phi_0;
        if (phi_0 >= 0) {
            m_phi_0 = pi - phi_0;
        } else {
            m_phi_0 = -phi_0 - pi;
        }
    } else {
        if (phi_0 <= pi) {
            m_phi_0 = phi_0;
        } else {
            m_phi_0 = phi_0 - pi2;
        }
    }
    
    psi_growth = part_growth * psi_end;
    psi_cutoff = part_cutoff * psi_end;
    
    mu_a_1 = log(r_end / r_0);
    m_a_1 = mu_a_1 - mu_a_2 + mu_a_3 - mu_a_4;
    m_a_2 = mu_a_2 - 3 * mu_a_3 + 6 * mu_a_4;
    m_a_3 = 2 * mu_a_3 - 10 * mu_a_4;
    m_a_4 = 5 * mu_a_4;
    
    bn = log(2);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double SpiralArm::GetValue(double x, double y) {
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

double SpiralArm::GetBrightness(double r, double psi) {
    if (psi <= 0 || psi >= psi_end)
        return 0;

    double r_spiral = GetRadius(psi);
    double I;
    I = I_0 * GetParallelBrightness(r_spiral, psi) * GetNormalBrightness(r_spiral, r - r_spiral, psi);

    return I;
}

double SpiralArm::GetParallelBrightness(double r_spiral, double psi) {
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

    return I_par;
}

double SpiralArm::GetNormalBrightness(double r_spiral, double rho, double psi) {
    double loc_w = abs(((w_i * r_spiral) + w_zp) / 2);
    return exp(-bn * pow(abs(rho) / loc_w, 2));
}

double SpiralArm::GetNearestCoordinates(double r, double psi) {
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

double SpiralArm::GetRadius(double psi) {
    double psi_norm = (psi / psi_end);
    return r_0 * exp(m_a_1 * psi_norm + m_a_2 * pow(psi_norm, 2) + m_a_3 * pow(psi_norm, 3) + m_a_4 * pow(psi_norm, 4));
}

/* END OF FILE: func_spiral.cpp ------------------------------------ */
