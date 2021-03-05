//
// Libraries
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//
// Structures
//

// Complex number
typedef struct {
    float re;    /* Real part */
    float im;    /* Imaginary part */
} complex_t;

//
// Main functions
//

// Sum
complex_t c_sum(complex_t z1, complex_t z2);

// Difference
complex_t c_diff(complex_t z1, complex_t z2);

// Product
complex_t c_product(complex_t z1, complex_t z2);

// Conjugate
complex_t c_conj(complex_t z);

// Module
float c_mod(complex_t z);

// Inner product
complex_t c_inner_product(complex_t z1, complex_t z2);

// Inner product between C3 vectors
complex_t c3_inner_product(complex_t *z1, complex_t *z2);

// Module
float c3_mod(complex_t *z);

// Module
float r3_mod(float *z);

// Inner product between R3 vectors
float r3_inner_product(float *r1, float *r2);

// Cross product between R3 vectors
float *r3_cross_product(float *r1, float *r2);

// Sum between R3 vectors
float *r3_sum(float *r1, float *r2);

// Difference between R3 vectors
float *r3_diff(float *r1, float *r2);

// Normalize a R3 vector
float *r3_normalize(float *r);
