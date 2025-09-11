/*   Class interface definition for func_spiral_0b.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * component with a broken-exponential profile.
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 *

 *
 *
 */

#include "function_object.h"

class SpiralArm0b : public FunctionObject {
    // the following static constant will be defined/initialized in the .cpp file
    static const char className[];

public:
    // Constructors:
    SpiralArm0b();

    // redefined method/member function:
    void Setup(double params[], int offsetIndex, double xc, double yc);

    double GetValue(double x, double y);
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName(string &classname) { classname = className; };


protected:
    double GetBrightness(double psi, double r);

    double GetNormalBrightness(double psi, double h);

    double GetNearestCoordinates(double r, double psi);

    double GetRadius(double psi);
private:
    double x0, y0, PA, ell, r_0, phi_0, r_end, phi_end, mu_a_2, mu_a_3, mu_a_4,
           I_0, part_growth, h_s, part_cutoff, width, w_asymm, n_out, n_in, gamma_out, gamma_in,
           q, cosPA, sinPA, is_clockwise, m_phi_0, psi_end, psi_growth, psi_cutoff,
           mu_a_1, m_a_1, m_a_2, m_a_3, m_a_4,
           n_out_inv, n_in_inv, bn, w_out, w_in;
};
