#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include "cosmical.h"
static char module_docstring[] = "This module allows Python to access Al Levine's fast and accurate cosmic ray generator written in C.";
static char cosmical_docstring[] = "Calculate a simulated cosmic ray image for a given exposure time.";

static PyObject *cosmical_cosmical(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"cosmical", cosmical_cosmical, METH_VARARGS, cosmical_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_cosmical(void)
{
    PyObject *m = Py_InitModule3("_cosmical", module_methods, module_docstring);
    if (m == NULL)
       return;

    /* Load 'numpy' functionality. */
    import_array();
}

//double cosmical(double exptm1, double exptm2, double crfl)


static PyObject *cosmical_cosmical(PyObject *self, PyObject *args)
{
    double crfl, exptm1, exptm2;
    long  NX, NY, idodif;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dddlll", &crfl, &exptm1, &exptm2, &NX, &NY, &idodif))
        return NULL;

    /* Call the external C function to compute the chi-squared. */
    double *unraveled = cosmical(crfl, exptm1, exptm2, NX, NY, idodif);

    npy_intp size[2];
    size[0] = NX;
    size[1] = NY;
    PyObject *ret = PyArray_SimpleNewFromData(2, size, NPY_DOUBLE, unraveled);
    PyArray_ENABLEFLAGS((PyArrayObject*) ret, NPY_ARRAY_OWNDATA);
    //Py_XDECREF(unraveled);

    /* Build the output tuple */
    //PyObject *ret = Py_BuildValue("d", py_array);
    return ret;
}
