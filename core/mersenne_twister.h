/** @file
 * \brief Pseudo-random-number-generator functions based on Mersenne Twister algorithm
 */

/* Header file for useful Mersenne Twister functions.
 * See mersenne_twister.c for full credits */

#ifndef _MERSENNE_TWISTER_H_
#define _MERSENNE_TWISTER_H_

/* Initialization routines: */
/*    initialize with a seed value */
void init_genrand( unsigned long s );
/*    initialize with an array */
void init_by_array( unsigned long init_key[], int key_length );

/* Integer rngs: */
/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32( );
/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31( );

/* Floating-point rngs: */
/* generates a random number on [0,1]-real-interval */
double genrand_real1( );
/* generates a random number on [0,1)-real-interval */
double genrand_real2( );
/* generates a random number on (0,1)-real-interval */
double genrand_real3( );
/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53( );


#endif /* _MERSENNE_TWISTER_H_ */
