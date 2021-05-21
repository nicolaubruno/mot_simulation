//
// Libraries
//

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "mot_sim.h"

//
// Functions
//

// Read simulation for a single atom
static PyObject* C_simulate_atom(PyObject *self, PyObject *args);

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void);

// Build a PyObject list of histograms
PyObject *build_freqs(histogram_t *pos);

// Build a PyObject list of a 3D-histogram
PyObject *build_3Dfreqs(histogram_3d_t pos);

//
// Structures
//

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