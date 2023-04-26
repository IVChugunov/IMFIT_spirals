/* FILE: func_flat-exp.cpp --------------------------------------------- */
/* VERSION 0.2
 *
 *   Function object class for an "flat/exponential" function, with constant
 * ellipticity and position angle (pure elliptical, not generalized).
 * The profile is an exponential for r > r_break, and constant intensity at
 * smaller radii.  This is a simplification of my broken-exponential function,
 * with the inner exponential slope set to infinity (to give a perfectly flat
 * inner profile).
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x and y.
 *      GetValue() then completes the calculation, using the actual value
 *      of x and y, and returns the result.
 *      So for an image, we expect the user to call Setup() once at
 *      the start, then loop through the pixels of the image, calling
 *      GetValue() to compute the function results for each pixel coordinate
 *      (x,y).
 *
 *   NOTE: Currently, we assume input PA is in *degrees* [and then we
 * convert it to radians] relative to +x axis.
 *
 *   MODIFICATION HISTORY:
 *     [v0.2]  31 Mar 2010: Tweaked to be an h1=infinity simplification of
 * the broken-exponential function, including variable smoothness of break.
 *     [v0.1]  30 Mar 2010: Created as modification of func_exp.cpp.
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_flat-exp.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_0", "h", "r_break", "alpha"};
const char  FUNCTION_NAME[] = "FlatExponential function";
//const char  SHORT_FUNCTION_NAME[] = "FlatExponential";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char FlatExponential::className[] = "FlatExponential";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

FlatExponential::FlatExponential( )
{
  string  paramName;
  
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

void FlatExponential::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_0 = params[2 + offsetIndex ];
  h = params[3 + offsetIndex ];
  r_b = params[4 + offsetIndex ];
  alpha = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  exponent = (1.0 / (alpha*h));
  double  S = pow( (1.0 + exp(-alpha*r_b)), exponent );
  I_0_times_S = I_0 * S;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for an exponential function at radius r,
// with the various parameters and derived values (n, b_n, r_e, etc.)
// pre-calculated by Setup().
// NOTE: We assume that r >= 0, since GetValue() ensures that.

double FlatExponential::CalculateIntensity( double r )
{
  double  I;
  
  // check for possible overflow in exponentiation if r >> r_b, and re-route around it:
  if ( alpha*(r - r_b) > 100.0 ) {
    // Outer-exponential approximation:
    I = I_0_times_S * exp(-r/h);
  } else {
    // no danger of overflow in exponentiation, so use fully correct calculation:
    I = I_0_times_S * pow( 1.0 + exp(alpha*(r - r_b)), -exponent );
  }
  return I;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double FlatExponential::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  
  nSubsamples = CalculateSubsamples(r);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-pixel
    double deltaSubpix = 1.0 / nSubsamples;
    double x_sub_start = x - 0.5 + 0.5*deltaSubpix;
    double y_sub_start = y - 0.5 + 0.5*deltaSubpix;
    double theSum = 0.0;
    for (int ii = 0; ii < nSubsamples; ii++) {
      double x_ii = x_sub_start + ii*deltaSubpix;
      for (int jj = 0; jj < nSubsamples; jj++) {
        double y_ii = y_sub_start + jj*deltaSubpix;
        x_diff = x_ii - x0;
        y_diff = y_ii - y0;
        xp = x_diff*cosPA + y_diff*sinPA;
        yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
        r = sqrt(xp*xp + yp_scaled*yp_scaled);
        theSum += CalculateIntensity(r);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// exponential function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int FlatExponential::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Based on standard exponential function subsampling, modified so that we
  // don't do subsampling if we're at r < r_break/2
  if ((doSubsampling) && (r > 0.5*r_b) && (r < 10.0)) {
    if ((h <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(2 * SUBSAMPLE_R / h));
    else {
      if (r <= 3.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}



/* END OF FILE: func_flat-exp.cpp -------------------------------------- */
