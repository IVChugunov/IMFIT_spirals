/*   Public interface for function which performs line-of-sight integration 
 * using Gnu Scientific Library's QAGS integration.
 */

#ifndef _GSL_INTEGRATOR_H_
#define _GSL_INTEGRATOR_H_

#include "gsl/gsl_integration.h"


double  Integrate( gsl_function F, double s1, double s2 );


#endif  // _GSL_INTEGRATOR_H_
