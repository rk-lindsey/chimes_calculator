#include <Python.h>




// Write the function that python will eventually call

static PyObject *method_fputs(PyObject *self, PyObject *args) 
{
   // Delcare argument types to recieve from python code

    char *str, *filename = NULL; 
    int bytes_copied = -1;

    /* Parse arguments */
    // "args" are of type PyObject
    // ss is the format specifier that specifies the data type of the arguments to parse
    // &str and &filename are pointers to local variable to which the parsed values will be assigned
    
    if(!PyArg_ParseTuple(args, "ss", &str, &filename)) 
    {
	return NULL;
    }
    
    // Standard C code you'd like to execute with python

    FILE *fp = fopen(filename, "w");
    bytes_copied = fputs(str, fp);
    fclose(fp);

    //PyLing_FromLong returns a PyLongObject which represents an integer in pytho 
    return PyLong_FromLong(bytes_copied);
}

// Write definitions of C module and the methods it contains... Include meta information about module
// tht will be used by thepython interpreter

// Tell python about your module (C code you'd like it to execute):
// fputs is the name the user would write to invoke function
// method_fputs is the name of the C function to invoke
// METH_VARARGS is a flag that tells the interpreter that the function will accepttwo arguments of type PyObject*
// Final string is a basic docstring
static PyMethodDef FputsMethods[] = {
    {"fputs", method_fputs, METH_VARARGS, "Python interface for fputs C library function"},
    {NULL, NULL, 0, NULL}
};

// PyModuleDef holds info on the C module itself... not an array of structures but rather a single structure
// used for module definition
// PyModuleDef_HEAD_INIT is usually never changed
// fputs is the name of the C extension module
// the string is a docstring
// -1 is the amount of memory needed to store the pogram state .. I don't know what that means
// FputsMethods is a reference to the method table, i.e. the array of PyMethodDef structus defined earlier
static struct PyModuleDef fputsmodule = {
    PyModuleDef_HEAD_INIT,
    "fputs",
    "Python interface for the fputs C library function",
    -1,
    FputsMethods
};

// Called by python the first time the C module is imported
// returns a new module object of type pyobject *
// for the argument, pass the address of the method structure previously defined, i.e. fputsmodule

PyMODINIT_FUNC PyInit_fputs(void) {
    return PyModule_Create(&fputsmodule);
}

