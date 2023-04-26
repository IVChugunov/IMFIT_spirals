/* FILE: integrator.cpp ------------------------------------------------ */
/* 
 * Code for performing a line-of-sight integration using GSL QAGS integration.
 * Everything is inside a single function, as local variables, to ensure thread
 * safety (e.g., for use with OpenMP).
 *
 * NOTE: Trial use of gsl_integration_qagi (integrating from -infty to +infty
 * sometimes worked, but sometimes failed (e.g., for 3D exponential disk when
 * inclination >~ 55 deg, though i = 90 worked).
 * Conclusion: safer to stick with gsl_integration_qags and do finite integration
 * to +/- large_number (s1 and s2).
 *
 * NOTE: Trial change of LIMIT_SIZE from 1000 to 10000, or RELATIVE_TOL from 1.0e-6 to
 * 1.0e-5, had no effect on integration failures for edges of moderately inclined 
 * Ferrers bar, so not much reason to change them.
 */
 
// Copyright 2011--2017 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include "integrator.h"

#define LIMIT_SIZE   1000
#define RELATIVE_TOL  1.0e-6

double  Integrate( gsl_function F, double s1, double s2 )
{
  double  result, error;
  int  status;
  gsl_integration_workspace * workspace;
  
  // allocate and free the workspace object here (referencing it with a local
  // variable) to ensure thread safety
  workspace = gsl_integration_workspace_alloc(LIMIT_SIZE);
  status = gsl_integration_qags(&F, s1, s2, 0, RELATIVE_TOL, LIMIT_SIZE, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  
  return result;
}


/* END OF FILE: integrator.cpp ----------------------------------------- */
