/** @file
 * \brief Utility functions for interpreting and printing results from fits.
 */

/* FILE: print_results.cpp ----------------------------------------- */

// Copyright 2010-2019 by Peter Erwin.
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


#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

#include "definitions.h"
#include "print_results.h"
#include "param_struct.h"
#include "statistics.h"
#include "utilities_pub.h"
#include "solver_results.h"
#ifndef NO_NLOPT
#include "nmsimplex_fit.h"
#include "nlopt_fit.h"
#endif

#define  FILE_OPEN_ERR_STRING "\n   Couldn't open file \"%s\"\n\n"



/* Local Functions: */

// Utility function for printing parameters to a file (including to stdout),
// with or without parameter errors (use NULL to indicate no errors) and
// with an optional prefix character for each line.
// Basically a wrapper around ModelObject::PrintModelParamsToStrings
void PrintParameters( FILE *filePtr, ModelObject *model, double *parameters, 
					  double *paramErrors, const char *prefix="" )
{
    vector<string> outputLines;

    model->PrintModelParamsToStrings(outputLines, parameters, paramErrors, prefix);
    for (string line: outputLines)
      fprintf(filePtr, "%s", line.c_str());
}

void GetSolverSummary( int status, int solverID, string& outputString );



// This is a function to print the results of a fit.  It's based on code from
// Craig Markwardt's testmpfit.c, but will also accomodate results from a fit
// done with other minimization algorithms, such as Nelder-Mead simplex or
// Differential Evolution (call with result=0 to indicate non-LM optimizer).
void PrintResults( double *params, ModelObject *model, int nFreeParameters, int fitStatus, 
					SolverResults& solverResults, bool recomputeStatistic )
{
  long  nValidPixels = model->GetNValidPixels();
  long  nDegreesFreedom = nValidPixels - nFreeParameters;
  int  whichStat, whichSolver;
  string  mpfitMessage;
  string  fitStatName, reducedStatName;
  double  fitStatistic, aic, bic;
  bool  printReduced;
  mp_result  *mpResult = nullptr;
  
  whichStat = model->WhichFitStatistic();
  whichSolver = solverResults.GetSolverType();

  if (recomputeStatistic)
    fitStatistic = model->GetFitStatistic(params);
  else
    fitStatistic = solverResults.GetBestfitStatisticValue();

  // Case of mpfit output (L-M minimization)
  if (whichSolver == MPFIT_SOLVER) {
    mpResult = solverResults.GetMPResults();
    InterpretMpfitResult(fitStatus, mpfitMessage);
    printf("\n*** mpfit status = %d -- %s\n", fitStatus, mpfitMessage.c_str());
    // Only print results of fit if valid fit was achieved
    if ((params == nullptr) || (mpResult == nullptr))
      return;
    if (whichStat == FITSTAT_CASH) {
      printf("  CASH STATISTIC = %f    (%ld DOF)\n", mpResult->bestnorm, nDegreesFreedom);
      printf("  INITIAL CASH STATISTIC = %f\n", mpResult->orignorm);
    }
    else if (whichStat == FITSTAT_POISSON_MLR) {
      printf("  POISSON-MLR STATISTIC = %f    (%ld DOF)\n", mpResult->bestnorm, nDegreesFreedom);
      printf("  INITIAL POISSON-MLR STATISTIC = %f\n", mpResult->orignorm);
    }
    else {
      printf("  CHI-SQUARE = %f    (%ld DOF)\n", mpResult->bestnorm, nDegreesFreedom);
      printf("  INITIAL CHI^2 = %f\n", mpResult->orignorm);
    }
    printf("        NPAR = %d\n", mpResult->npar);
    printf("       NFREE = %d\n", mpResult->nfree);
    printf("     NPEGGED = %d\n", mpResult->npegged);
    printf("     NITER = %d\n", mpResult->niter);
    printf("      NFEV = %d\n", mpResult->nfev);
    printf("\n");
    aic = AIC_corrected(mpResult->bestnorm, nFreeParameters, nValidPixels, 1);
    bic = BIC(mpResult->bestnorm, nFreeParameters, nValidPixels, 1);
    if (whichStat == FITSTAT_CHISQUARE)
      printf("Reduced Chi^2 = %f\n", mpResult->bestnorm / nDegreesFreedom);
    if (whichStat == FITSTAT_POISSON_MLR)
      printf("Reduced Chi^2 equivalent = %f\n", mpResult->bestnorm / nDegreesFreedom);
    printf("AIC = %f, BIC = %f\n", aic, bic);
    
    double *paramErrs = (double *)calloc(model->GetNParams(), sizeof(double));
    solverResults.GetErrors(paramErrs);
    PrintParameters(stdout, model, params, paramErrs);
    printf("\n");

    free(paramErrs);
  }
  else {
    // Only print results of fit if fitStatus >= 1
    if (fitStatus < 1)
      return;
    fitStatName = "CHI-SQUARE";
    reducedStatName = "Reduced Chi^2";
    printReduced = true;
    if (whichStat == FITSTAT_CASH) {
      fitStatName = "CASH STATISTIC";
      printReduced = false;
    }
    else if (whichStat == FITSTAT_POISSON_MLR) {
      fitStatName = "POISSON-MLR STATISTIC";
      reducedStatName = "Reduced Chi^2 equivalent";
    }
    printf("\n  %s = %f\n", fitStatName.c_str(), fitStatistic);
    printf("  NFEV = %d\n\n", solverResults.GetNFunctionEvals());
    if (printReduced)
      printf("%s = %f\n", reducedStatName.c_str(), fitStatistic / nDegreesFreedom);
    aic = AIC_corrected(fitStatistic, nFreeParameters, nValidPixels, 1);
    bic = BIC(fitStatistic, nFreeParameters, nValidPixels, 1);
    printf("AIC = %f, BIC = %f\n", aic, bic);
    PrintParameters(stdout, model, params, NULL);
    printf("\n");
  }
}

/// Simple function for computing and printing the fit statistic given the input
/// parameter vector. Meant for use in imfit_main.cpp with the "--fitstat-only"
/// (aka "chisquare-only") option
void PrintFitStatistic( double *params, ModelObject *model, int nFreeParameters )
{
  long  nValidPixels = model->GetNValidPixels();
  long  nDegreesFreedom = nValidPixels - nFreeParameters;
  int  whichStat = model->WhichFitStatistic();
  double  fitStatistic = model->GetFitStatistic(params);
  double  aic = AIC_corrected(fitStatistic, nFreeParameters, nValidPixels, 1);
  double  bic = BIC(fitStatistic, nFreeParameters, nValidPixels, 1);

  if (whichStat == FITSTAT_CASH)
    printf("  CASH STATISTIC = %f\n", fitStatistic);
  else if (whichStat == FITSTAT_POISSON_MLR) {
    printf("  POISSON-MLR STATISTIC = %f\n", fitStatistic);
    printf("\nReduced Chi^2 equivalent = %f\n", fitStatistic / nDegreesFreedom);
  }
  else {
    printf("  CHI-SQUARE = %f    (%ld DOF)\n", fitStatistic, nDegreesFreedom);
    printf("\nReduced Chi^2 = %f\n", fitStatistic / nDegreesFreedom);
  }
  printf("AIC = %f, BIC = %f\n\n", aic, bic);

}



void GetSolverSummary( int status, int solverID, string& outputString )
{
  string  tempString;
  
  outputString = "Algorithm: ";
  switch (solverID) {
    case MPFIT_SOLVER:
      outputString += PrintToString("Levenberg-Marquardt: status = %d -- ", status);
      InterpretMpfitResult(status, tempString);
      outputString += tempString;
      break;
#ifndef NO_NLOPT
    case NMSIMPLEX_SOLVER:
      GetInterpretation_NM(status, tempString);
      outputString += tempString;
      break;
    case GENERIC_NLOPT_SOLVER:
      GetInterpretation_NLOpt(status, tempString);
      outputString += tempString;
      break;
#endif
    case DIFF_EVOLN_SOLVER:
      outputString += PrintToString("Differential Evolution: status = %d -- ", status);
      if (status == 1)
        outputString += "SUCCESS: Convergence in fit-statistic value";
      else  // assuming (status == 5)
        outputString += "Maximum generation number reached without convergence";
      break;
  }
}


/// Saves best-fit parameters (and summary of fit statistics) to a file.
void SaveParameters( double *params, ModelObject *model, string& outputFilename, 
					vector<string>& outputHeader, int nFreeParameters, int whichSolver, 
					int fitStatus, SolverResults& solverResults )
{
  FILE  *file_ptr;
  string  statName, algorithmSummary;
  double  *parameterErrs = NULL;
  
  if ((file_ptr = fopen(outputFilename.c_str(), "w")) == NULL) {
    fprintf(stderr, FILE_OPEN_ERR_STRING, outputFilename.c_str());
    exit(-1);
  }

  // Get minimization (solver output) info
  GetSolverSummary(fitStatus, whichSolver, algorithmSummary);
  
  // Get fit-results info
  long  nValidPixels = model->GetNValidPixels();
  long  nDegreesFreedom = nValidPixels - nFreeParameters;
  double  fitStatistic = model->GetFitStatistic(params);
  double  aic = AIC_corrected(fitStatistic, nFreeParameters, nValidPixels, 1);
  double  bic = BIC(fitStatistic, nFreeParameters, nValidPixels, 1);
  int  whichStat = model->WhichFitStatistic();
  if (whichStat == FITSTAT_CASH) {
    statName = "Cash statistic";
  }
  else if (whichStat == FITSTAT_POISSON_MLR) {
    statName = "Poisson-MLR statistic";
  }
  else {
    whichStat = model->WhichFitStatistic(true);
    switch (whichStat) {
      case FITSTAT_CHISQUARE_USER:
        statName = "chi-squared (user-supplied errors)";
        break;
      case FITSTAT_CHISQUARE_MODEL:
        statName = "chi-squared (model-based errors)";
        break;
      default:
        statName = "chi-squared (data-based errors)";
        break;
    }
  }
  
  for (int i = 0; i < (int)outputHeader.size(); i++)
    fprintf(file_ptr, "%s\n", outputHeader[i].c_str());
  fprintf(file_ptr, "\n# Results of fit:\n");
  fprintf(file_ptr, "#   %s\n", algorithmSummary.c_str());
  fprintf(file_ptr, "#   Fit statistic: %s\n", statName.c_str());
  fprintf(file_ptr, "#   Best-fit value: %f\n", fitStatistic);
  if (whichStat == FITSTAT_CASH) {
    fprintf(file_ptr, "#   Reduced value: none\n");
  }
  else {
    fprintf(file_ptr, "#   Reduced value: %f\n", fitStatistic / nDegreesFreedom);
  }
  fprintf(file_ptr, "#   AIC: %f\n", aic);
  fprintf(file_ptr, "#   BIC: %f\n", bic);

  if (solverResults.ErrorsPresent()) {
    parameterErrs = (double *)calloc(model->GetNParams(), sizeof(double));
    solverResults.GetErrors(parameterErrs);
    PrintParameters(file_ptr, model, params, parameterErrs);
    free(parameterErrs);
  }
  else
    PrintParameters(file_ptr, model, params, NULL);

  fclose(file_ptr);
}



// Same as SaveParameters, but requires a previously opened file pointer; also allows
// optional prefix string for each line (default declaration in header file = "")
void SaveParameters2( FILE *file_ptr, double *params, ModelObject *model, 
                    vector<string>& outputHeader, const char *prefix )
{
  for (int i = 0; i < (int)outputHeader.size(); i++)
    fprintf(file_ptr, "%s\n", outputHeader[i].c_str());
  fprintf(file_ptr, "%s\n", prefix);
  PrintParameters(file_ptr, model, params, NULL, prefix);
}


/* END OF FILE: print_results.cpp ---------------------------------- */
