//
// Libraries
//

#include <Python.h>
#include "mot_sim.h"

//
// Implementation of C Extension for Python
//

// Read the parameters from CSV files and generate the results
static PyObject* run(PyObject *self, PyObject *args){
    return Py_BuildValue("i", C_run());
}

// List methods
static PyMethodDef mot_sim_methods[] = {
    {"run", run, METH_NOARGS, "Read the parameters from CSV files and generate the results"},
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