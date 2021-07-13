//
// Libraries
//

#include "mot_sim.h"
#include <time.h>

int main(){
    int i, j;
    results_t res;
    
    res = simulate_atom("../parameters/", 2, time(0));
    printf("trapped_atom = %d\n", res.trapped_atom);
    printf("time [tau] = %f\n", res.time);


    return 0;
}