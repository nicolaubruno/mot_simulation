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
// Complex Space
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

// Sum between C3 vectors
complex_t *c3_sum(complex_t *z1, complex_t *z2);

// Difference between C3 vectors
complex_t *c3_diff(complex_t *z1,complex_t *z2);

// Product between a complex number and a real vector
complex_t *c3_scalar_product(complex_t a, complex_t *v);

// Inner product between C3 vectors
complex_t c3_inner_product(complex_t *z1, complex_t *z2);

// Module
float c3_mod(complex_t *z);

// View complex vector
int c3_view(complex_t *z, char *name);

// Apply complex operator
complex_t *c3_apply_operator(complex_t **A, complex_t *v);

// Convert a R3 vector to a C3 vector
complex_t *r3_to_c3(float *v);

//
// Real space
//

// Module
float r3_mod(float *z);

// Inner product between R3 vectors
float r3_inner_product(float *r1, float *r2);

// Cross product between R3 vectors
float *r3_cross_product(float *r1, float *r2);

// Product between a real number and a real vector
float *r3_scalar_product(float a, float *v);

// Sum between R3 vectors
float *r3_sum(float *r1, float *r2);

// Difference between R3 vectors
float *r3_diff(float *r1, float *r2);

// Normalize a R3 vector
float *r3_normalize(float *r);

// Apply real operator
float *r3_apply_operator(float **A, float *v);

// View real vector
int r3_view(float *z, char *name);
