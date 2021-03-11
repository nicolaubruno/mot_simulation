//
// Libraries
//

#include "mot_sim.h"

int main(){
    int i, j, k;
    results_t res;

    for(i = 0; i < 10; i++){
        res = simulate_atom("../results/1615476589_test/loop0/parameters/");
        printf("%d\n", res.num_iters);
    }

    return 0;
}