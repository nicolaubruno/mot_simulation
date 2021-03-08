//
// Libraries
//

#include <Python.h>
#include "mot_sim.h"

//
// Implementation of C Extension for Python
//

// Read simulation for a single atom
static PyObject* C_simulate_atom(PyObject *self, PyObject *args){
    //
    // Variables
    //
    char *dir_code;

    // Parse arguments
    if(!PyArg_ParseTuple(args, "s", &dir_code)) return NULL;

    return Py_BuildValue("i", simulate_atom(dir_code));
}

// List methods
static PyMethodDef mot_sim_methods[] = {
    {"simulate_atom", C_simulate_atom, METH_VARARGS, "Read the parameters from CSV files and generate the results"},
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

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void) {
    return PyModule_Create(&mot_sim_module);
}