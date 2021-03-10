//
// Libraries
//

#include <Python.h>
#include "mot_sim.h"

//
// Prototypes
//

// Read simulation for a single atom
static PyObject* C_simulate_atom(PyObject *self, PyObject *args);

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void);

// Convert a int array in C to PyObject list
PyObject *build_list(int *arr, int size);

//
// Structures
//

// List methods
static PyMethodDef mot_sim_methods[] = {
    {"simulate_atom", C_simulate_atom, METH_NOARGS, "Read the parameters from CSV files and generate the results"},
    {NULL, NULL, 0, NULL}
};

// Define module
static struct PyModuleDef mot_sim_module = {
    PyModuleDef_HEAD_INIT,
    "mot_sim",
    "C Extension with optimized functions related to the MOTSim",
    -1,
    mot_sim_methods
};

//
// Implementation
//

// Read simulation for a single atom
static PyObject* C_simulate_atom(PyObject *self, PyObject *args){
    //
    // Results
    //

    results_t res = simulate_atom();
    PyObject *x_bins, *y_bins, *z_bins, *ret;

    // Build PyLists
    x_bins = build_list(res.pos_hist[0].freqs, res.pos_hist[0].num_bins);
    y_bins = build_list(res.pos_hist[1].freqs, res.pos_hist[1].num_bins);
    z_bins = build_list(res.pos_hist[2].freqs, res.pos_hist[2].num_bins);

    // Return
    ret = Py_BuildValue(
        "OOOid", 
        x_bins, 
        y_bins,
        z_bins,
        res.num_iters,
        res.time
    );

    Py_DECREF(x_bins);
    Py_DECREF(y_bins);
    Py_DECREF(z_bins);

    return ret;
}

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void){
    return PyModule_Create(&mot_sim_module);
}

// Convert a int array in C to PyObject list
PyObject *build_list(int *arr, int size){
    //
    // Variables
    //

    int i;
    PyObject *list;

    //
    // Build list
    //

    list = PyList_New(size);

    for(i = 0; i < size; i++){
        PyList_SetItem(list, i, Py_BuildValue("i", arr[i]));
    }

    // Return
    return list;
}