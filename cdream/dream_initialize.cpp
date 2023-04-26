#include <cmath>
#include <iostream>
#include <vector>
#include <stdio.h>
using namespace std;

#include "dream.h"

void dream_initialize( const dream_pars* p, rng::RngStream* rng, 
						Array2DView<double>& state, ArrayView<double>& lik) 
{
  if (p->verboseLevel > 0)
    printf("Latin hypercube sampling...\n");
  
  Array2D<int> samples(p->nvar, p->numChains);
  samples.set_all(0);

  if (p->verboseLevel > 0)
    printf("  shuffling (%d, %d)...", p->nvar, p->numChains);
  for (int j = 0; j < p->nvar; ++j) {
    for (int i = 0; i < p->numChains; ++i) 
      samples(j,i) = i;
    rng->shuffle(samples(j), p->numChains);
  }
  if (p->verboseLevel > 0)
    printf("done.\n");

  double rand;
  double randPos;
  double intervalSize;

  for (int i = 0; i < p->numChains; ++i) {
    for (int j = 0; j < p->nvar; ++j) {
      if (p->varLock[j]) 
        randPos = p->varLo[j];
      else {
        rng->uniform(1, &rand);
        intervalSize = (p->varHi[j] - p->varLo[j])/p->numChains;
        randPos = p->varLo[j] + samples(j,i)*(p->varHi[j] - p->varLo[j])/p->numChains
                    + intervalSize*rand;
      }
      if (randPos < p->varLo[j]) 
        randPos = p->varLo[j];
      else if (randPos > p->varHi[j]) 
        randPos = p->varHi[j];
      state(i,j) = randPos;
      if (p->verboseLevel > 1) 
        printf("(%d,%d) = %f\n", i, j, state(i,j));
    }
  }

  int do_calc = 1;
//  const double* vars = NULL;

  for (int i = 0; i < p->numChains; ++i) {
    do_calc = 1;
    for (int j = 0; j < p->nvar; ++j) {
      if (state(i,j) < p->varLo[j] || state(i,j) > p->varHi[j]) {
        do_calc = 0;
        break;
      }
    }
    if (do_calc) {
      lik[i] = p->fun(i, -1, state.col_pt(i), p->extraData, false);
      if (p->verboseLevel > 0)
        printf("Chain %d: likelihood = %.2f\n", i, lik[i]);
//        cout << "Chain " << i << " likelihood = " << flush;
//      if (p->verboseLevel > 0)
//        cout << lik[i] << endl;
    } else 
      lik[i] = -INFINITY;
  }

  // replace chains with infinite likelihood with random samples
  for (int i = 0; i < p->numChains; ++i) {
    while (lik[i] == -INFINITY || lik[i] != lik[i]) { // TODO: ?????
      printf("Likelihood of chain = INF. Resampling intial parameters...\n");
      // choose random parameters
      for (int j = 0; j < p->nvar; ++j) {
        if (p->varLock[j]) 
          randPos = p->varLo[j];
        else {
          rng->uniform(1, &rand);
          randPos = p->varLo[j] + (p->varHi[j]-p->varLo[j])*rand;
        }
        if (randPos < p->varLo[j]) 
          randPos = p->varLo[j];
        else if (randPos > p->varHi[j]) 
          randPos = p->varHi[j];
        state(i,j) = randPos;
      }
      // if (fixedRatio >= 0.0) state(0,i,1) = state(0,i,2)*(1./fixedRatio-1.);
      lik[i] = p->fun(i, -1, state.col_pt(i), p->extraData, false);
      printf("New likelihood = %f\n", lik[i]);
    }
  }
}


