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

#include <random>

#include <stdio.h>
#include "dcdio.h"

#include "smc.h"

using namespace std;

extern void get_distances();

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

float energy(float *x_array, float *y_array, float *z_array, energy_parameters p) {

    /*
        method to calculate energy.  
        returns the sum of the hard-sphere and long-range energy
    */ 
    int i, j ;
    float overlap = 1E99 ;
    float u_long_range = 0.0 ;
    float u_long_range_1, u_long_range_2 ;

    float xi, yi, zi, xf, yf, zf, dx, dy, dz, dx2, dy2, dz2, r ;

    for (i = 0 ; i < p.natoms - 1 ; i ++)
    {
        xi = x_array[i];
        yi = y_array[i];
        zi = z_array[i];
    
        for (j = i + 1 ; j < p.natoms ; j++)
        {
            xf = x_array[j];
            yf = y_array[j];
            zf = z_array[j];
   
            dx = xf - xi ; dy = yf - yi ; dz = zf - zi ;  
            dx2 = dx * dx ; dy2 = dy * dy ; dz2 = dz * dz ;
            r = sqrt(dx2 + dy2 + dz2) ;

            if (r < ( p.sigma_1 + p.sigma_2)) {
                return overlap ;
            }
 
            u_long_range_1 = -p.epsilon_ab_a * pow(( p.sigma_ab / p.rab_a ),2.0) * exp(-(r/p.rab_a)) ;
            u_long_range_2 =  p.epsilon_ab_r * pow(( p.sigma_ab / p.rab_r ),2.0) * exp(-(r/p.rab_r)) ;

            u_long_range += u_long_range_1 + u_long_range_2 ; 

        } // end of j-loop

    } // end of i-loop


    return u_long_range ;
}


int surface_move(float *x_array, float * y_array, float *z_array, int i, energy_parameters parameters) {

    float x, y, z, r, dx, dy, dz, tx, ty, tz ;
    float max_disp, norm ;
    float u_long_range, boltz, delta_energy, ran ; 

    std::random_device rd ;
    std::random_device rd2 ;
    std::mt19937 mt(rd());
    std::mt19937 mt2(rd2());
    std::uniform_real_distribution<float> dist(-1.0,1.0);
    std::uniform_real_distribution<float> dist2(0.0,1.0);

    bool accepted = false ;

    x = x_array[i] ; 
    y = y_array[i] ; 
    z = z_array[i] ; 

    tx = x ; ty = y ; tz = z ;
    //std::cout << "c : x = " << x << std::endl ;
    
    r = sqrt(x*x+y*y+z*z) ;

    max_disp = 0.05 ;
    //max_disp = 1.0 ;
    
    dx = max_disp * dist(mt);
    dy = max_disp * dist(mt);
    dz = max_disp * dist(mt);

    x += dx ; y += dy ; z += dz ;
    norm = sqrt(x*x + y*y + z*z) ;
    x *= r/norm ;
    y *= r/norm ;
    z *= r/norm ;

    x_array[i] = x ;
    y_array[i] = y ;
    z_array[i] = z ;

    u_long_range = energy(x_array, y_array, z_array, parameters) ;

    if (u_long_range < parameters.energy) {
        accepted = true ;
    } else {
        delta_energy = u_long_range - parameters.energy ;
        boltz = exp(-parameters.beta * delta_energy) ;
        ran = dist2(mt2) ;
        if (ran < boltz) {
            accepted = true ;
        } 
    }
    if (accepted) {
        parameters.energy = u_long_range ;
        return 1 ;
    } else {
        x_array[i] = tx ;
        y_array[i] = ty ;
        z_array[i] = tz ;
        return 0 ;
    }
}

PyObject *smc_parallel(PyObject *self, PyObject *args){
	PyObject *array = NULL ;
    PyObject *pList = NULL ;

    //std::ostringstream sstream;
    //std::string remark = "#hello";
    //const std::string filename = "dum.txt"; 

    //std::ofstream outfile(filename.c_str()) ;
    //outfile << remark << std::endl;
    int i, j, k ; 
    int number_of_steps, natoms ;
    double temperature, sigma_1, sigma_2, epsilon_ab_a, epsilon_ab_r, rab_a, rab_r, sigma_ab, beta, contrast_1, contrast_2 ;

    const char * dcdfile_name ;

    if (!PyArg_ParseTuple(args, "OOiiddddddddddds", &array, &pList, &number_of_steps, &natoms, &temperature, &sigma_1, &sigma_2, &epsilon_ab_a, &epsilon_ab_r, &rab_a, &rab_r, &sigma_ab, &beta, &contrast_1, &contrast_2, &dcdfile_name))
        return NULL;

    std::cout << "c: number_of_steps = " << number_of_steps << std::endl ; 
    std::cout << "c: natoms = " <<  natoms<< std::endl ; 
    std::cout << "c: temperature = " <<  temperature << std::endl ; 
    std::cout << "c: sigma_1 = " <<  sigma_1 << std::endl ; 
    std::cout << "c: sigma_2 = " <<  sigma_2 << std::endl ; 
    std::cout << "c: epsilon_ab_a = " <<  epsilon_ab_a << std::endl ; 
    std::cout << "c: epsilon_ab_r = " <<  epsilon_ab_r << std::endl ; 
    std::cout << "c: rab_a = " <<  rab_a << std::endl ; 
    std::cout << "c: rab_r = " <<  rab_r << std::endl ; 
    std::cout << "c: sigma_ab = " <<  sigma_ab << std::endl ; 
    std::cout << "c: beta = " <<  beta << std::endl ; 
    std::cout << "c: contrast_1 = " <<  contrast_1 << std::endl ; 
    std::cout << "c: contrast_2 = " <<  contrast_2 << std::endl ; 
    std::cout << "c: dcdfile_name = " <<  dcdfile_name << std::endl ; 


    // put energy inputs into an instance of a struct
    //
    //typedef struct energy_parameters parameters ;
    struct energy_parameters parameters ;
    parameters.natoms = natoms ;
    parameters.temperature = temperature ;
    parameters.sigma_1 = sigma_1 ;
    parameters.sigma_2 = sigma_2 ;
    parameters.epsilon_ab_a = epsilon_ab_a ;
    parameters.epsilon_ab_r = epsilon_ab_r ;
    parameters.rab_a = rab_a ;
    parameters.rab_r = rab_r ;
    parameters.sigma_ab = sigma_ab ;
    parameters.beta = beta ;
    parameters.energy = 1E99 ;

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
    for(i = 0 ; i < natoms ; i++){
        //printf("c : %i\t %f \n", i,c_array[0][i][0]);
        x_array[i] = (float)c_array[0][i][0] ;
        y_array[i] = (float)c_array[0][i][1] ;
        z_array[i] = (float)c_array[0][i][2] ;
    }

    printf("c : %f \n", c_array[0][0][0]);
    std::cout << "c : x_array[0] = " << x_array[0] << std::endl ;
    std::cout << "c : y_array[0] = " << y_array[0] << std::endl ;
    std::cout << "c : z_array[0] = " << z_array[0] << std::endl ;

    int *c_int_array;
    int typenum2 = NPY_INT;
    PyArray_Descr *descr2;
    descr2 = PyArray_DescrFromType(typenum2);
    npy_intp dims2[1];
    if (PyArray_AsCArray(&pList, (void *)&c_int_array, dims2, 1, descr2) < 0) {
        PyErr_SetString(PyExc_TypeError, "error converting to c array");
        return NULL;
    }

    std::cout << "c: pList[0] = " << c_int_array[0] << std::endl ;
    std::cout << "c: pList[1] = " << c_int_array[1] << std::endl ;

    int number_accepted = 0 ;

    // Open DCD file for writing 
    //
    //
    char * filename = (char*)dcdfile_name ;

    FILE * filepointer  ;
    filepointer = open_dcd_write(filename) ;

    int istart = 0 ;
    int nsavc = 1 ;
    int nset = 1 ;
    double delta = 1.0 ;
    int headerresult ;
    int stepresult ;
    headerresult = write_dcdheader(filepointer, filename, natoms, nset, istart, nsavc, delta) ;

    std::cout << "write dcd header result = " << headerresult << std::endl ;

    for(i = 0 ; i < number_of_steps ; i++) {
     //   std::cout << "c : x_array[0] = " << x_array[0] << std::endl ;
        std::cout << i << std::flush ;   
        for(j = 0 ; j < natoms ; j++){ 
            number_accepted += surface_move(x_array, y_array, z_array, j, parameters) ;
        }

    //    std::cout << "c : x_array[0] = " << x_array[0] << std::endl ;
            // stepresult = dcdio.write_dcdstep(filepointer, tx, ty, tz, step)¬
        stepresult = write_dcdstep(filepointer, natoms, x_array, y_array, z_array, i) ;

    } // end of i loop

    stepresult = close_dcd_write(filepointer) ;

    std::cout << "\n\nI am ready to return\n\n" << std::endl ;

    //get_distances(c_array, nframes, natoms, hist, nbins, bin_width) ;

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

