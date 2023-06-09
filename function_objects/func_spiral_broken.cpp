
/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_spiral_broken.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int N_PARAMS = 24;
const char PARAM_LABELS[][30] = {"m0", "m1", "m2", "m3", "fi_break", "l0", "l1", "l2", "l3", "PA", "inc", "fi0", "r0", "fi_max", "n", "is_clockwise",
                                 "I0", "r_e", "width_i", "outer_r_e", "inner_r_e", "fi_of_max", "outer_n", "inner_n"};
const char FUNCTION_NAME[] = "Spiral branch broken function";
const double pi = 3.14159265358979323846;
const double pi2 = 6.28318530718;

const char SpiralBranchBroken::className[] = "SpiralBranchBroken";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

SpiralBranchBroken::SpiralBranchBroken() {
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

void SpiralBranchBroken::Setup(double params[], int offsetIndex, double xc, double yc) {
    x0 = xc;
    y0 = yc;
    m0 = params[0 + offsetIndex];
    m1 = params[1 + offsetIndex];
    m2 = params[2 + offsetIndex];
    m3 = params[3 + offsetIndex];
    fi_break = params[4 + offsetIndex] / 180 * pi;
    l0 = params[5 + offsetIndex];
    l1 = params[6 + offsetIndex];
    l2 = params[7 + offsetIndex];
    l3 = params[8 + offsetIndex];
    PA = params[9 + offsetIndex] / 180 * pi;
    inc = params[10 + offsetIndex] / 180 * pi;
    fi0 = params[11 + offsetIndex] / 180 * pi;
    r0 = params[12 + offsetIndex];
    fi_max = params[13 + offsetIndex] / 180 * pi;
    n = 1 / (params[14 + offsetIndex] / 180 * pi);
    is_clockwise = params[15 + offsetIndex];
    max_bright = params[16 + offsetIndex];
    bright_decrease = params[17 + offsetIndex];
    width_increase = params[18 + offsetIndex];
    outer_width = params[19 + offsetIndex];
    inner_width = params[20 + offsetIndex];
    fi_of_max = params[21 + offsetIndex] / 180 * pi;
    outer_n = params[22 + offsetIndex];
    inner_n = params[23 + offsetIndex];

    outer_bn = Calculate_bn(outer_n);
    inner_bn = Calculate_bn(inner_n);

    mp = bright_decrease;
    kp = abs(GetRadius(fi_of_max) - r0) / mp;

    maximum = max_bright / (pow(mp * kp, kp) * exp(-kp));
    if (kp <= 1e-4)
        maximum = max_bright;

    if (is_clockwise > 0) {
        f_fi0 = pi - fi0;
    } else {
        if (fi0 <= pi) {
            f_fi0 = fi0;
        } else {
            f_fi0 = fi0 - 2 * pi;
        }
    }
    oinvn = 1.0 / outer_n;
    iinvn = 1.0 / inner_n;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double SpiralBranchBroken::GetValue(double x, double y) {
    double xc = x - x0;
    double yc = y - y0;

    double xp = xc * cos(PA) + yc * sin(PA);
    double yp = -xc * sin(PA) + yc * cos(PA);

    double x_d0 = xp;
    double y_d0 = yp / cos(inc);

    double r = sqrt(x_d0 * x_d0 + y_d0 * y_d0);

    if (is_clockwise > 0)
        x_d0 = -x_d0;

    double fi = atan2l(y_d0, x_d0) - f_fi0;

    if (fi < 0)
        fi += pi * 2;
    if (fi >= pi * 2)
        fi -= pi * 2;

    double fiMax = fi_max + fi_break + 1 / n;
    double minFi = GetNearestCoordinates(r, fi);
    double maxFi;
    if (minFi == 0)
        maxFi = fi;
    else
        maxFi = minFi + 2 * pi;

    while (minFi > fiMax)
        minFi -= 2 * pi;

    double bright = 0;
    bright += GetBrightness(minFi, r);
    bright += GetBrightness(maxFi, r);
    return bright;
}

/* ----------------------------- OTHER FUNCTIONS -------------------------------- */

double SpiralBranchBroken::GetBrightness(double fi, double r) {
    if (fi <= 0 || fi >= fi_max + fi_break + 1 / n)
        return 0;

    double spiral_r = GetRadius(fi);
    double bright = maximum * pow(abs(spiral_r - r0), kp) * exp(-abs(spiral_r - r0) / mp) * GetNormalBrightness(fi, r - spiral_r);

    if (fi > fi_max + fi_break) {
        return bright * (1 - n * (fi - (fi_max + fi_break)));
    } else
        return bright;
}


double SpiralBranchBroken::GetNormalBrightness(double fi, double h) {
    if (fi == 0)
        return 0;
    double width = h < 0 ? inner_width : outer_width;
    double invn = h < 0 ? iinvn : oinvn;
    double bn = h < 0 ? inner_bn : outer_bn;
    return exp(-bn * pow(abs(h) / sqrt(width * width + pow(width_increase * fi, 3)), invn) );
}

double SpiralBranchBroken::GetRadius(double fi) {
    double x = fi / pi2;
    double x_break = fi_break / pi2;
    double r_break = r0 * exp((m0 + m1 * x_break + m2 * pow(x_break, 2) + m3 * pow(x_break, 3)) * fi_break);
    double fi_on2 = fi - fi_break;
    double x_on2 = fi_on2 / pi2;
    if (fi > fi_break) {
    	return r_break * exp((l0 + l1 * x_on2 + l2 * pow(x_on2, 2) + l3 * pow(x_on2, 3)) * fi_on2);
    } else
    	return r0 * exp((m0 + m1 * x + m2 * pow(x, 2) + m3 * pow(x, 3)) * fi);
}

double SpiralBranchBroken::GetNearestCoordinates(double r, double fi) {
    double r2 = GetRadius(fi);
    double fi1 = 0;
    double fi2 = fi;
    int i = 0;
    while (r2 < r && i++ < 15) {
        fi1 = fi2;
        fi2 = fi2 + pi * 2;
        r2 = GetRadius(fi2);
    }

    return fi1;
}



/* END OF FILE: func_spiral_broken.cpp ------------------------------------ */
