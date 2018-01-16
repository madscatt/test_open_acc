#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"
#include "oacc_smc.h"
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


/** Convert a c++ 2D vector into a numpy array
 *
 * @param const vector< vector<T> >& vec : 2D vector data
 * @return PyArrayObject* array : converted numpy array
 *
 * Transforms an arbitrary 2D C++ vector into a numpy array. Throws in case of
 * unregular shape. The array may contain empty columns or something else, as
 * long as it's shape is square.
 *
 * Warning this routine makes a copy of the memory!
 */
template<typename T> static PyArrayObject* vector_to_nparray(const vector< vector<T> >& vec, int type_num = PyArray_INT){

   // rows not empty
   if( !vec.empty() ){

      // column not empty
      if( !vec[0].empty() ){

        size_t nRows = vec.size();
        size_t nCols = vec[0].size();
        npy_intp dims[2] = {nRows, nCols};
        PyArrayObject* vec_array = (PyArrayObject *) PyArray_SimpleNew(2, dims, type_num);

        T *vec_array_pointer = (T*) PyArray_DATA(vec_array);

        // copy vector line by line ... maybe could be done at one
        for (size_t iRow=0; iRow < vec.size(); ++iRow){

          if( vec[iRow].size() != nCols){
             Py_DECREF(vec_array); // delete
             throw(string("Can not convert vector<vector<T>> to np.array, since c++ matrix shape is not uniform."));
          }

          copy(vec[iRow].begin(),vec[iRow].end(),vec_array_pointer+iRow*nCols);
        }

        return vec_array;

     // Empty columns
     } else {
        npy_intp dims[2] = {vec.size(), 0};
        return (PyArrayObject*) PyArray_ZEROS(2, dims, PyArray_INT, 0);
     }

   // no data at all
   } else {
      npy_intp dims[2] = {0, 0};
      return (PyArrayObject*) PyArray_ZEROS(2, dims, PyArray_INT, 0);
   }

}


PyObject *smc_parallel(PyObject *self, PyObject *args){
	PyObject *array = NULL ;
    PyObject *pylist, *item ;

    //std::ostringstream sstream;
    //std::string remark = "#hello";
    //const std::string filename = "dum.txt"; 

    //std::ofstream outfile(filename.c_str()) ;
    //outfile << remark << std::endl;
     
    int nframes, natoms, nbins ;
    double bin_width ;

    if (!PyArg_ParseTuple(args, "Oiiid", &array, &nframes, &natoms, &nbins, &bin_width))
        return NULL;

    std::cout << "c: nbins = " << nbins << std::endl ; 
    std::cout << "c: bin_width = " << bin_width << std::endl ; 

    double ***c_array;
    //unsigned long long int npairs ;
    //std::vector<double> dist(npairs, 0.0) ;
    //std::vector<int> hist(nbins, 0) ;
    std::vector<std::vector<int> > hist(nframes, std::vector<int>(nbins,       0));
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

    get_distances(c_array, nframes, natoms, hist, nbins, bin_width) ;

    PyObject *np_vec_2D = (PyObject*) vector_to_nparray(hist) ;

    printf("leaving smc_parallel\n") ;

    return np_vec_2D ;
}

static PyMethodDef methods[] = {
	{ "smc_parallel", smc_parallel, METH_VARARGS },
	{ NULL}
} ;

extern "C" void initsmc_parallel(){
	PyObject *m ;
	m = Py_InitModule("smc_parallel", methods);
	import_array();
}

