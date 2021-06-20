//
// Libraries
//

#include "mot_sim.h"
#include <time.h>

int main(){
    int i, j;
    results_t res;

    for(i = 0; i < 1; i++){
        res = simulate_atom("../parameters/", 2, time(0));
        
        //printf("\nSim %d\n", i+1);

        //printf("x = [");
        //for(j = 0; j < res.pos_hist[0].num_bins; j++)
            //printf(" %d ", res.pos_hist[0].freqs[j]);
        //printf("]\n");

        //printf("\n");

        //printf("y = [");
        //for(j = 0; j < res.pos_hist[1].num_bins; j++)
        //    printf(" %d ", res.pos_hist[1].freqs[j]);
        //printf("]\n");

        //printf("\n");

        //printf("z = [");
        //for(j = 0; j < res.pos_hist[2].num_bins; j++)
        //    printf(" %d ", res.pos_hist[2].freqs[j]);
        //printf("]\n\n");
    }


    return 0;
}