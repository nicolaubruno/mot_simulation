//
// Libraries
//

#include "Py_wrapper.h"

// Read simulation for a single atom
static PyObject* C_simulate_atom(PyObject *self, PyObject *args){
    //
    // Variables
    int i, j;
    int only_marginals;
    long time;
    PyObject *ret;
    results_t res;
    char *params_path;
    PyObject *pos_freqs, *vel_freqs;

    //
    // Get parameters path
    if(!PyArg_ParseTuple(args, "sil", &params_path, &only_marginals, &time)){
        return NULL;
    }
    
    //
    // Results
    res = simulate_atom(params_path, only_marginals, time);
    //
    // Position and velocity frequencies
    //--

    // Only marginal distributions
    if(only_marginals){
        // Build values
        pos_freqs = build_freqs(res.pos_hist);
        vel_freqs = build_freqs(res.vel_hist);

        //
        // Release memory
        for(i = 0; i < 3; i++){
            free(res.pos_hist[i].freqs);
            free(res.vel_hist[i].freqs);
        }

        free(res.pos_hist);
        free(res.vel_hist);
    // 3D-Distributions
    } else {
        // Build values
        pos_freqs = build_3Dfreqs(res.pos_3Dhist);
        vel_freqs = build_3Dfreqs(res.vel_3Dhist);

        //
        // Release memory
        //--
        for(i = 0; i < res.pos_3Dhist.num_bins[0]; i++){
            for(j = 0; j < res.pos_3Dhist.num_bins[1]; j++){
                free(res.pos_3Dhist.freqs[i][j]);
            }

            free(res.pos_3Dhist.freqs[i]);
        }
        
        free(res.pos_3Dhist.freqs);
        free(res.pos_3Dhist.num_bins);
        free(res.pos_3Dhist.bins_size);
        free(res.pos_3Dhist.coord0);


        for(i = 0; i < res.vel_3Dhist.num_bins[0]; i++){
            for(j = 0; j < res.vel_3Dhist.num_bins[1]; j++){
                free(res.vel_3Dhist.freqs[i][j]);
            }

            free(res.vel_3Dhist.freqs[i]);
        }
        
        free(res.vel_3Dhist.freqs);
        free(res.vel_3Dhist.num_bins);
        free(res.vel_3Dhist.bins_size);
        free(res.vel_3Dhist.coord0);
        //--
    }
    //--

    // Return
    ret = Py_BuildValue("OOi", pos_freqs, vel_freqs, time);

    return ret;
}

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void){
    return PyModule_Create(&mot_sim_module);
}

// Convert the results in a PyObject list
PyObject *build_3Dfreqs(histogram_3d_t hist){
    //
    // Variables
    //

    int i, j, k, dim_x, dim_y, dim_z;
    PyObject *list;

    //
    // Build list
    //

    dim_x = hist.num_bins[0];
    dim_y = hist.num_bins[1];
    dim_z = hist.num_bins[2];

    list = PyList_New(dim_x*dim_y*dim_z);

    for(i = 0; i < dim_x; i++){
        for(j = 0; j < dim_y; j++){
            for(k = 0; k < dim_z; k++){
                PyList_SetItem(list, (dim_y*dim_z)*i + dim_z*j + k, Py_BuildValue("i", hist.freqs[i][j][k]));
            }
        }
    }

    // Return
    return list;
}

// Convert the results in a PyObject list
PyObject *build_freqs(histogram_t *hist){
    //
    // Variables
    //

    int i, j, dim;
    PyObject *list, *item, *item2;

    //
    // Build list
    //

    list = PyList_New(3);

    for(i = 0; i < 3; i++){
        dim = hist[i].num_bins;
        item = PyList_New(dim);

        for(j = 0; j < dim; j++){
            item2 = Py_BuildValue("i", hist[i].freqs[j]);
            PyList_SetItem(item, j, item2);
        }
        
        PyList_SetItem(list, i, item);
    }

    // Return
    return list;
}