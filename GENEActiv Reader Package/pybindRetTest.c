#include <Python.h>

/*
 * Implements an example function.
 */
PyDoc_STRVAR(pybindRetTest_example_doc, "example(obj, number)\
\
Example function");

PyObject *pybindRetTest_example(PyObject *self, PyObject *args, PyObject *kwargs) {
    /* Shared references that do not need Py_DECREF before returning. */
    PyObject *obj = NULL;
    int number = 0;

    /* Parse positional and keyword arguments */
    static char* keywords[] = { "obj", "number", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", keywords, &obj, &number)) {
        return NULL;
    }

    /* Function implementation starts here */

    if (number < 0) {
        PyErr_SetObject(PyExc_ValueError, obj);
        return NULL;    /* return NULL indicates error */
    }

    Py_RETURN_NONE;
}

/*
 * List of functions to add to pybindRetTest in exec_pybindRetTest().
 */
static PyMethodDef pybindRetTest_functions[] = {
    { "example", (PyCFunction)pybindRetTest_example, METH_VARARGS | METH_KEYWORDS, pybindRetTest_example_doc },
    { NULL, NULL, 0, NULL } /* marks end of array */
};

/*
 * Initialize pybindRetTest. May be called multiple times, so avoid
 * using static state.
 */
int exec_pybindRetTest(PyObject *module) {
    PyModule_AddFunctions(module, pybindRetTest_functions);

    PyModule_AddStringConstant(module, "__author__", "maczak");
    PyModule_AddStringConstant(module, "__version__", "1.0.0");
    PyModule_AddIntConstant(module, "year", 2024);

    return 0; /* success */
}

/*
 * Documentation for pybindRetTest.
 */
PyDoc_STRVAR(pybindRetTest_doc, "The pybindRetTest module");


static PyModuleDef_Slot pybindRetTest_slots[] = {
    { Py_mod_exec, exec_pybindRetTest },
    { 0, NULL }
};

static PyModuleDef pybindRetTest_def = {
    PyModuleDef_HEAD_INIT,
    "pybindRetTest",
    pybindRetTest_doc,
    0,              /* m_size */
    NULL,           /* m_methods */
    pybindRetTest_slots,
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL,           /* m_free */
};

PyMODINIT_FUNC PyInit_pybindRetTest() {
    return PyModuleDef_Init(&pybindRetTest_def);
}
