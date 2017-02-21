#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"
//#include "oacc_pr.h"
#include <vector> 

#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


using namespace std;

extern void get_distances();

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/


PyObject *cm_parallel(PyObject *self, PyObject *args){
	PyObject *array_1 = NULL ;
	PyObject *array_2 = NULL ;
    PyObject *pylist, *item ;

    std::ostringstream sstream;
    std::string remark = "#hello";
    const std::string filename = "dum.txt"; 

    std::ofstream outfile(filename.c_str()) ;
    outfile << remark << std::endl;
     
    int nframes, natoms_1, natoms_2, nbins ;
    double bin_width ;
    double cutoff ;

    if (!PyArg_ParseTuple(args, "OOiiiidd", &array_1, &array_2, &nframes, &natoms_1, &natoms_2, &nbins, &bin_width, &cutoff))
        return NULL;

    std::cout << "c: nbins = " << nbins << std::endl ; 
    std::cout << "c: bin_width = " << bin_width << std::endl ; 

    double ***c_array_1;
    double ***c_array_2;

    unsigned long long int npairs ;
    //std::vector<double> dist(npairs, 0.0) ;
    std::vector<int> hist(nbins, 0) ;


    //Create C arrays from numpy objects:
    int typenum = NPY_DOUBLE;
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(typenum);
    npy_intp dims[3];

    if (PyArray_AsCArray(&array_1, (void ***)&c_array_1, dims, 3, descr) < 0) {
        PyErr_SetString(PyExc_TypeError, "error converting to c array");
        return NULL;
    }
    
    if (PyArray_AsCArray(&array_2, (void ***)&c_array_2, dims, 3, descr) < 0) {
        PyErr_SetString(PyExc_TypeError, "error converting to c array");
        return NULL;
    }

    printf("c : nframes = %d\n", nframes);
    printf("c : natoms_1 = %d\n", natoms_1);
    printf("c : natoms_2 = %d\n", natoms_2);
    printf("c : %f \n", c_array_1[0][0][0]);
    printf("c : %f \n", c_array_2[0][0][0]);
    printf("c : cutoff = %f \n", cutoff);
/*
    get_distances(c_array, nframes, natoms, hist, nbins, bin_width) ;

    //pylist = PyList_New(npairs) ;
    pylist = PyList_New(nbins) ;
    if (pylist != NULL){
        for (int i=0 ; i<nbins ; i++) {
            sstream << hist[i] << endl;
            //std::string value = sstream.str();
            item = PyInt_FromLong(hist[i]);
            //item = PyFloat_FromDouble(hist[i]);
            //PyList_SET_ITEM(pylist, i, item);
            PyList_SetItem(pylist, i, item);
        }
    }
    outfile << sstream.str() ;

    //free dist ;

    return pylist ;
 */
    return Py_None ;

}

static PyMethodDef methods[] = {
	{ "cm_parallel", cm_parallel, METH_VARARGS },
	{ NULL}
} ;

extern "C" void initcm_parallel(){
	PyObject *m ;
	m = Py_InitModule("cm_parallel", methods);
	import_array();
}

