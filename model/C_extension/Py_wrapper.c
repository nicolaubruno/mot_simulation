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

    //
    // Get parameters path
    if(!PyArg_ParseTuple(args, "sil", &params_path, &only_marginals, &time)){
        return NULL;
    }
    
    //
    // Results
    res = simulate_atom(params_path, only_marginals, time);

    if(only_marginals){
        // Build values
        ret = Py_BuildValue("O", build_pos_freqs(res));

        //
        // Release memory
        for(i = 0; i < 3; i++){
            free(res.pos_hist[i].freqs);
        }

        free(res.pos_hist);
    } else {
        // Build values
        ret = Py_BuildValue("O", build_pos_3Dfreqs(res));

        //
        // Release memory
        //

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
    }

    // Release memory
    //free(params_path);

    return ret;
}

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void){
    return PyModule_Create(&mot_sim_module);
}

// Convert the results in a PyObject list
PyObject *build_pos_3Dfreqs(results_t res){
    //
    // Variables
    //

    int i, j, k, dim_x, dim_y, dim_z;
    PyObject *list;

    //
    // Build list
    //

    dim_x = res.pos_3Dhist.num_bins[0];
    dim_y = res.pos_3Dhist.num_bins[1];
    dim_z = res.pos_3Dhist.num_bins[2];

    list = PyList_New(dim_x*dim_y*dim_z);

    for(i = 0; i < dim_x; i++){
        for(j = 0; j < dim_y; j++){
            for(k = 0; k < dim_z; k++){
                PyList_SetItem(list, (dim_y*dim_z)*i + dim_z*j + k, Py_BuildValue("i", res.pos_3Dhist.freqs[i][j][k]));
            }
        }
    }

    // Return
    return list;
}

// Convert the results in a PyObject list
PyObject *build_pos_freqs(results_t res){
    //
    // Variables
    //

    int i, j, dim;
    PyObject *list, *item;

    //
    // Build list
    //

    list = PyList_New(3);

    for(i = 0; i < 3; i++){
        dim = res.pos_hist[i].num_bins;
        item = PyList_New(dim);

        for(j = 0; j < dim; j++)
            PyList_SetItem(item, j, Py_BuildValue("i", res.pos_hist[i].freqs[j]));
        
        PyList_SetItem(list, i, item);
    }

    // Return
    return list;
}