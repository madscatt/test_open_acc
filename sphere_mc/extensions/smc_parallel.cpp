#include "Python.h"
#include "numpy/arrayobject.h"
#include <iostream>

#ifndef SMC_H
#define SMC_H

#include "smc.h"

#endif

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

PyObject *smc_parallel(PyObject *self, PyObject *args){
	PyObject *array = NULL ;
    PyObject *pList = NULL ;

    int i ;
    int natoms, number_of_steps, ncell_1d, ncell ;
    double temperature, sigma_11, sigma_22, sigma_12 ;
    double epsilon_a_11, epsilon_a_22, epsilon_a_12 ;
    double epsilon_r_11, epsilon_r_22, epsilon_r_12 ;
    double r_a_11, r_a_22, r_a_12, r_r_11, r_r_22, r_r_12 ; 
    double beta, contrast_1, contrast_2, r_cutoff, max_displacement;
    float cell_length, delta ;
    const char * dcdfile_name ;

    if (!PyArg_ParseTuple(args, "OOffiiiiddddddddddddddddddddds", &array, &pList, &cell_length, &delta, &ncell_1d, &ncell, &number_of_steps, &natoms, &temperature, &sigma_11, &sigma_22, &sigma_12,  &epsilon_a_11, &epsilon_a_22, &epsilon_a_12, &epsilon_r_11, &epsilon_r_22, &epsilon_r_12, &r_a_11, &r_a_22, &r_a_12, &r_r_11, &r_r_22, &r_r_12, &beta, &contrast_1, &contrast_2, &r_cutoff, &max_displacement, &dcdfile_name))
        return NULL;

    std::cout << "c: number_of_steps = " << number_of_steps << std::endl ; 
    std::cout << "c: natoms = " <<  natoms<< std::endl ; 
    std::cout << "c: temperature = " <<  temperature << std::endl ; 
    std::cout << "c: sigma_11 = " <<  sigma_11 << std::endl ; 
    std::cout << "c: sigma_22 = " <<  sigma_22 << std::endl ; 
    std::cout << "c: sigma_12 = " <<  sigma_12 << std::endl ; 
    std::cout << "c: epsilon_a_11 = " <<  epsilon_a_11 << std::endl ; 
    std::cout << "c: epsilon_a_22 = " <<  epsilon_a_22 << std::endl ; 
    std::cout << "c: epsilon_a_12 = " <<  epsilon_a_12 << std::endl ; 
    std::cout << "c: epsilon_r_11 = " <<  epsilon_r_11 << std::endl ; 
    std::cout << "c: epsilon_r_22 = " <<  epsilon_r_22 << std::endl ; 
    std::cout << "c: epsilon_r_12 = " <<  epsilon_r_12 << std::endl ; 
    std::cout << "c: r_a_11 = " <<  r_a_11 << std::endl ; 
    std::cout << "c: r_a_22 = " <<  r_a_22 << std::endl ; 
    std::cout << "c: r_a_12 = " <<  r_a_12 << std::endl ; 
    std::cout << "c: r_r_11 = " <<  r_r_11 << std::endl ; 
    std::cout << "c: r_r_22 = " <<  r_r_22 << std::endl ; 
    std::cout << "c: r_r_12 = " <<  r_r_12 << std::endl ; 
    std::cout << "c: beta = " <<  beta << std::endl ; 
    std::cout << "c: contrast_1 = " <<  contrast_1 << std::endl ; 
    std::cout << "c: contrast_2 = " <<  contrast_2 << std::endl ; 
    std::cout << "c: r_cutoff = " <<  r_cutoff << std::endl ; 
    std::cout << "c: max_displacement = " <<  max_displacement << std::endl ; 
    std::cout << "c: dcdfile_name = " <<  dcdfile_name << std::endl ; 

    // put energy inputs into an instance of a struct
    //
    struct energy_parameters parameters ;
    parameters.natoms = natoms ;
    parameters.temperature = (float)temperature ;
    parameters.sigma_11 = (float)sigma_11 ;
    parameters.sigma_22 = (float)sigma_11 ;
    parameters.sigma_12 = (float)sigma_12 ;
    parameters.epsilon_a_11 = (float)epsilon_a_11;
    parameters.epsilon_a_22 = (float)epsilon_a_22;
    parameters.epsilon_a_12 = (float)epsilon_a_12;
    parameters.epsilon_r_11 = (float)epsilon_r_11 ;
    parameters.epsilon_r_22 = (float)epsilon_r_22 ;
    parameters.epsilon_r_12 = (float)epsilon_r_12 ;
    parameters.r_a_11 = (float)r_a_11 ;
    parameters.r_a_22 = (float)r_a_22 ;
    parameters.r_a_12 = (float)r_a_12 ;
    parameters.r_r_11 = (float)r_r_11 ;
    parameters.r_r_22 = (float)r_r_22 ;
    parameters.r_r_12 = (float)r_r_12 ;
    parameters.beta = (float)beta ;
    parameters.energy = (float)1E99 ;
    parameters.r_cutoff = (float)r_cutoff ;
    parameters.max_displacement = (float)max_displacement ;
    parameters.contrast_1 = (float)contrast_1 ;
    parameters.contrast_2 = (float)contrast_2 ;

    // put system inputs into an instance of a struct
    //
    // since some values are const, assignment is done in order
    struct system_parameters system_parameters  = {\
        number_of_steps,\
        cell_length,\
        delta,\
        ncell_1d,\
        ncell,\
        dcdfile_name\
    } ;

    double ***c_array;

    float x_array[natoms] ;
    float y_array[natoms] ;
    float z_array[natoms] ;

    //Create C arrays from numpy objects:
    int typenum = NPY_DOUBLE;
    PyArray_Descr *descr;
    descr = PyArray_DescrFromType(typenum);
    npy_intp dims[3];
    if (PyArray_AsCArray(&array, (void ***)&c_array, dims, 3, descr) < 0) {
        PyErr_SetString(PyExc_TypeError, "error converting to c array");
        return NULL;
    }

    // store coordinates in sepearte 1D arrays
    //
    for(i = 0 ; i < natoms ; ++i){
        //printf("c : %i\t %f \n", i,c_array[0][i][0]);
        x_array[i] = (float)c_array[0][i][0] ;
        y_array[i] = (float)c_array[0][i][1] ;
        z_array[i] = (float)c_array[0][i][2] ;
    }

    printf("c : %f \n", c_array[0][0][0]);
    std::cout << "c : x_array[0] = " << x_array[0] << std::endl ;
    std::cout << "c : y_array[0] = " << y_array[0] << std::endl ;
    std::cout << "c : z_array[0] = " << z_array[0] << std::endl ;

    int *atom_id;
    int typenum2 = NPY_INT;
    PyArray_Descr *descr2;
    descr2 = PyArray_DescrFromType(typenum2);
    npy_intp dims2[1];
    if (PyArray_AsCArray(&pList, (void *)&atom_id, dims2, 1, descr2) < 0) {
        PyErr_SetString(PyExc_TypeError, "error converting to c array");
        return NULL;
    }

    std::cout << "c: pList[0] = " << atom_id[0] << std::endl ;
    std::cout << "c: pList[1] = " << atom_id[1] << std::endl ;

    smc_core(x_array, y_array, z_array, atom_id, parameters, system_parameters) ;

    Py_RETURN_NONE ;

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

