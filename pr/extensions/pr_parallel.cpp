#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"
#include "oacc_pr.h"
#include <vector> 
using namespace std;

extern void get_distances();

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/


PyObject *pr_parallel(PyObject *self, PyObject *args){
	PyObject *array = NULL ;
    PyObject *pylist, *item ;
    
    int nframes, natoms ;

    if (!PyArg_ParseTuple(args, "Oii", &array, &nframes, &natoms))
        return NULL;

    double ***c_array;
    unsigned long long int npairs ;
    npairs = (natoms * (natoms - 1))/2 ;
    double dist[(natoms * (natoms -1))/2] ;
    //std:vector<double> dist(npairs, 0.0) ;

    //Create C arrays from numpy objects:
    int typenum = NPY_DOUBLE;
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(typenum);
    npy_intp dims[3];
    if (PyArray_AsCArray(&array, (void ***)&c_array, dims, 3, descr) < 0) {
        PyErr_SetString(PyExc_TypeError, "error converting to c array");
        return NULL;
    }

    printf("c : nframes = %d\n", nframes);
    printf("c : natoms = %d\n", natoms);
    printf("c : %f \n", c_array[0][0][0]);

    get_distances(c_array, nframes, natoms, dist) ;

    //for(int i=0 ; i < 30 ; i++){
     //   printf("dist[%i] = %f\n", i,dist[i]/double(nframes)) ;
   // } 

    pylist = PyList_New(npairs) ;
    if (pylist != NULL){
        for (unsigned long long int i=0 ; i<npairs ; i++) {
            item = PyFloat_FromDouble(dist[i]/double(nframes));
            PyList_SET_ITEM(pylist, i, item);
        }
    }

    //free dist ;

    return pylist ;
   // return Py_None ;

}

static PyMethodDef methods[] = {
	{ "pr_parallel", pr_parallel, METH_VARARGS },
	{ NULL}
} ;

extern "C" void initpr_parallel(){
	PyObject *m ;
	m = Py_InitModule("pr_parallel", methods);
	import_array();
}

