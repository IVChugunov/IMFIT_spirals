/*   Class interface definition for func_spiral.cpp
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

class SpiralBranch : public FunctionObject {
    // the following static constant will be defined/initialized in the .cpp file
    static const char className[];

public:
    // Constructors:
    SpiralBranch();

    // redefined method/member function:
    void Setup(double params[], int offsetIndex, double xc, double yc);

    double GetValue(double x, double y);
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName(string &classname) { classname = className; };


protected:
    double GetBrightness(double fi, double r);

    double GetNormalBrightness(double fi, double h);

    double GetNearestCoordinates(double r, double fi);

    double GetRadius(double fi);
private:
    double x0, y0, is_clockwise, PA, inc, r0, fi0, max_bright, bright_decrease, width_increase,
            outer_width, inner_width, m0, m1, m2, m3, fi_max, n, fi_of_max, kp, mp, maximum, cp,
            outer_n, inner_n, outer_bn, inner_bn, oinvn, iinvn;
    double f_fi0;
};
