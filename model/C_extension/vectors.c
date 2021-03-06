//
// Libraries
//

#include "vectors.h"

//
// Complex space
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

// Sum between C3 vectors
complex_t *c3_sum(complex_t *z1, complex_t *z2){
    //
    // Variables
    //

    int i;
    static complex_t res[3];

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = c_sum(z1[i], z2[i]);

    return res;
}

// Difference between C3 vectors
complex_t *c3_diff(complex_t *z1, complex_t *z2){
    //
    // Variables
    //

    int i;
    static complex_t res[3];

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = c_diff(z1[i], z2[i]);

    return res;
}

// Product between a complex number and a real vector
complex_t *c3_scalar_product(complex_t a, complex_t *v){
    // Variables
    int i;
    static complex_t *w;

    w = (complex_t*) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++)
        w[i] = c_product(a, v[i]);

    return w;
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

// View complex vector
int c3_view(complex_t *z, char *name){
    int i;

    printf("%s = [", name);
    for(i = 0; i < 2; i++) printf("(%f, %f), ", z[i].re, z[i].im);
    printf("(%f, %f)]\n", z[i].re, z[i].im);
}

// Apply complex operator
complex_t *c3_apply_operator(complex_t **A, complex_t *v){
    int i,  j;
    static complex_t *a;

    a = (complex_t*) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            a[i] = c_sum(a[i], c_product(A[i][j], v[j]));
        }
    }

    return a;
}

// Convert a R3 vector to a C3 vector
complex_t *r3_to_c3(float *v){
    int i;
    static complex_t *res;

    res = (complex_t *) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++){
        res[i].re = v[i];
        res[i].im = 0;
    }

    return res;
}

//
// Real space
//

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
    for(i = 0; i < 3; i++) res += r1[i] * r2[i];

    return res;
}

// Cross product between R3 vectors
float *r3_cross_product(float *r1, float *r2){
    //
    // Variables
    //

    int i;
    static float *res;

    res = (float*) calloc(3, sizeof(float));

    // Compute dot product
    res[0] = r1[1]*r2[2] - r1[2]*r2[1];
    res[1] = r1[2]*r2[0] - r1[0]*r2[2];
    res[2] = r1[0]*r2[1] - r1[1]*r2[0];

    return res;
}

// Product between a real number and a real vector
float *r3_scalar_product(float a, float *v){
    // Variables
    int i;
    static float w[3];

    for(i = 0; i < 3; i++)
        w[i] = a * v[i];

    return w;
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

// Apply real operator
float *r3_apply_operator(float **A, float *v){
    int i,  j;
    static float *a;

    a = (float*) calloc(3, sizeof(float));

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            a[i] += A[i][j]*v[j];
        }
    }

    return a;
}


// View real vector
int r3_view(float *z, char *name){
    int i;

    printf("%s = [", name);
    for(i = 0; i < 2; i++) printf("%f, ", z[i]);
    printf("%f]\n", z[i]);
}