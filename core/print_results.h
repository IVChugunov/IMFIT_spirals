/** @file
 * \brief Public interface for functions which print the results of a fit 
 *        to stdout or to file.
 *
 *    Utility functions for interpreting and printing results from fits
 *
 */

#ifndef _PRINT_RESULTS_H_
#define _PRINT_RESULTS_H_

#include <string>
#include "model_object.h"
#include "param_struct.h"
#include "solver_results.h"


/// Code for printing the results of a fit to stdout.
void PrintResults( double *params, ModelObject *model, int nFreeParameters, int fitStatus, 
					SolverResults& solverResults, bool recomputeStatistic=false );

void PrintFitStatistic( double *params, ModelObject *model, int nFreeParameters );

void SaveParameters( double *params, ModelObject *model, string& outputFilename, 
					vector<string>& outputHeader, int nFreeParameters, int whichSolver, 
					int fitStatus, SolverResults& solverResults );

void SaveParameters2( FILE *file_ptr, double *params, ModelObject *model, 
                    	vector<string>& outputHeader, const char *prefix );


#endif /* _PRINT_RESULTS_H_ */
