//
// Libraries
//

#include "vectors.h"

//
// Functions
//

// Sum
complex_t c_sum(complex_t z1, complex_t z2){
    complex_t z;

    z.re = z1.re + z2.re;
    z.im = z1.im + z2.im;

    return z;
}

// Difference
complex_t c_diff(complex_t z1, complex_t z2){
    complex_t z;

    z.re = -z2.re;
    z.im = -z2.im;

    return c_sum(z1, z);
}

// Product
complex_t c_product(complex_t z1, complex_t z2){
    complex_t z;

    z.re = z1.re * z2.re - z1.im * z2.im;
    z.im = z1.re * z2.im + z1.im * z2.re;

    return z; 
}

// Conjugate
complex_t c_conj(complex_t z){
    z.im = -z.im;
    return z;
}

// Inner product
complex_t c_inner_product(complex_t z1, complex_t z2){
    return c_product(z1, c_conj(z2));
}

// Module
float c_mod(complex_t z){
    return sqrt(c_inner_product(z, z).re);
}

// Inner product between C3 vectors
complex_t c3_inner_product(complex_t *z1, complex_t *z2){
    int i;
    complex_t z;

    z.re = 0;
    z.im = 0;

    for(i = 0; i < 3; i++){
        z = c_sum(z, c_inner_product(z1[i], z2[i]));
    }

    return z;
}

// Module
float c3_mod(complex_t *z){
    return sqrt(c3_inner_product(z, z).re);
}

// Module
float r3_mod(float *z){
    return sqrt(r3_inner_product(z, z));
}


// Inner product between R3 vectors
float r3_inner_product(float *r1, float *r2){
    //
    // Variables
    //

    int i;
    float res = 0;

    // Compute dot product
    for(i = 0; i < 3; i++) res += (*(r1+i)) * (*(r2+i));

    return res;
}

// Cross product between R3 vectors
float *r3_cross_product(float *r1, float *r2){
    //
    // Variables
    //

    int i;
    static float res[3];

    // Compute dot product
    res[0] = r1[1]*r2[2] - r1[2]*r2[1];
    res[1] = r1[2]*r2[0] - r1[0]*r2[2];
    res[2] = r1[0]*r2[1] - r1[1]*r2[0];

    return res;
}

// Sum between R3 vectors
float *r3_sum(float *r1, float *r2){
    //
    // Variables
    //

    int i;
    static float res[3];

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = r1[i] + r2[i];

    return res;
}

// Subtraction between R3 vectors
float *r3_diff(float *r1, float *r2){
    //
    // Variables
    //

    int i;
    static float res[3];

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = r1[i] - r2[i];

    return res;
}

// Normalize a R3 vector
float *r3_normalize(float *r){
    //
    // Variables
    //

    static float *new_r;
    float mod_r;
    int i;

    new_r = (float*) calloc(3, sizeof(float));

    //
    // Normalization
    //

    mod_r = sqrt(r3_inner_product(r, r));
    if(mod_r == 0.0){
        printf("Division by zero in normalize_vector() function\n");
        exit(0);
    }
    for(i = 0; i < 3; i++) new_r[i] = (r[i] / mod_r);

    return new_r;
}