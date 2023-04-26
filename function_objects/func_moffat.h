/*   Class interface definition for func_moffat.cpp
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * Moffat function.
 *
 * PARAMETERS:
 * PA = params[0 + offsetIndex];
 * ell = params[1 + offsetIndex];
 * I_0 = params[2 + offsetIndex ];
 * fwhm = params[3 + offsetIndex ];
 * beta = params[4 + offsetIndex ];
 *
 *
 */


// CLASS Moffat:

#include "function_object.h"



/// Class for image function with elliptical isophotes and %Moffat profile
class Moffat : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Moffat( );
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
    double  x0, y0, PA, ell, I_0, fwhm, beta;   // parameters
    double  alpha;
    double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};
