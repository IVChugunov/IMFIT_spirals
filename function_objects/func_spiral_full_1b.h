/*   Class interface definition for func_spiral_full_1b.cpp
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

class SpiralArmFull1b : public FunctionObject {
    // the following static constant will be defined/initialized in the .cpp file
    static const char className[];

public:
    // Constructors:
    SpiralArmFull1b();

    // redefined method/member function:
    void Setup(double params[], int offsetIndex, double xc, double yc);

    double GetValue(double x, double y);
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName(string &classname) { classname = className; };


protected:
    double GetBrightness(double r, double psi);

    double GetParallelBrightness(double r_spiral, double psi);

    double GetNormalBrightness(double r, double rho, double psi);

    double GetNearestCoordinates(double r, double psi);

    double GetRadius(double psi);
private:
    double x0, y0, PA, ell, r_0, phi_0, r_break_1, phi_break_1, r_end, phi_end, mu_a_2, mu_a_3, mu_b_2, mu_b_3,
           I_0, part_growth, ih_s, part_cutoff, w_zp, w_i, S_0, S_end, n,
           q, cosPA, sinPA, is_clockwise, m_phi_0, psi_break_1, psi_end, psi_growth, psi_cutoff,
           mu_a_1, mu_b_1, m_a_1, m_a_2, m_a_3, m_b_1, m_b_2, m_b_3,
           n_inv, bn, S_zp, S_i,
           A_d_1, phi_d_1, w_d_1, A_d_2, phi_d_2, w_d_2,
           psi_d_1, psi_d_2;
};
