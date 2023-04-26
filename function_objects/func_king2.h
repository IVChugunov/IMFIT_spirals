/*   Class interface definition for func_king.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for a modified
 * King model (e.g., Elson 1999; Peng+2010).
 *   This is the variant version with concentration as a free parameter
 * (instead of r_t)
 *
 * UNFINISHED!
 *
 * PARAMETERS:
 * x0 = xc;   -- center of component (pixels, x)
 * y0 = yc;   -- center of component (pixels, y)
 * PA = params[0 + offsetIndex];   -- PA of component, rel. to +x axis
 * ell = params[1 + offsetIndex];  -- ellipticity
 * I_0 = params[2 + offsetIndex ]; -- central intensity (ADU/pix)
 * r_c = params[3 + offsetIndex ];   -- core radius (pixels)
 * c   = params[4 + offsetIndex ];   -- concentration
 * alpha = params[5 + offsetIndex ];   -- power-law index
 *
 *
 */


// CLASS ModifiedKing2:

#include "function_object.h"
#include <string>
using namespace std;


/// Class for image function with elliptical isophotes and modified King profile (using concentration as free parameter instead of tidal radius)
class ModifiedKing2 : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    ModifiedKing2( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
   // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };

  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, I_0, r_c, c, alpha;   // parameters
    double  q, PA_rad, cosPA, sinPA;   // other useful, geometry-related quantities
    double  r_t, I_1, one_over_alpha, one_over_rc, constantTerm;   // other useful, profile-related quantities
};

