#include <Python.h>
#include "sim_atom.h"

static PyObject * factorial(PyObject *self, PyObject *args){
    int n;

    if(!PyArg_ParseTuple(args, "i", &n))
        return NULL;

    return Py_BuildValue("i", C_factorial(n));
}

static PyMethodDef sim_atom_methods[] = {
    {"factorial", factorial, METH_VARARGS, "Calculate the factorial of an integer number"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef sim_atom_module = {
    PyModuleDef_HEAD_INIT,
    "sim_atom",
    "Simulation of a single atom",
    -1,
    sim_atom_methods
};

PyMODINIT_FUNC PyInit_sim_atom(void) {
    return PyModule_Create(&sim_atom_module);
}