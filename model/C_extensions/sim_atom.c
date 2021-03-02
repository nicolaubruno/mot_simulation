#include "sim_atom.h"

int C_factorial(int n){
    if (n > 1) {
        return n * C_factorial(n - 1);
    } else {
        return 1;
    }
}