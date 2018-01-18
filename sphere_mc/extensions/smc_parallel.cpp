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

using namespace std;

extern void get_distances();

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/


int surface_move(float *x_array, float * y_array, float *z_array, int i) {

    float x, y, z, r, dx, dy, dz ;;
    float max_disp, norm ;

    std::random_device rd ;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(-1.0,1.0);

    bool accepted = false ;

    x = x_array[i] ; 
    y = y_array[i] ; 
    z = z_array[i] ; 

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

    //std::cout << dx << std::endl ;
    //std::cout << dy << std::endl ;
    //std::cout << dz << std::endl << std::endl ;
/*
    u_long_range = energy(mol, p)

    if u_long_range:
        if u_long_range < p.energy:
            accepted = true ;
        else:
            delta_energy = u_long_range - p.energy
            boltz = math.exp(-p.beta * delta_energy)
            ran = random.random()
            if ran < boltz:
                accepted = true ;

        if accepted:
            p.energy = u_long_range
            mol.coor()[0][i][0] = x
            mol.coor()[0][i][1] = y
            mol.coor()[0][i][2] = z
*/
    //std::cout << "c : x = " << x << std::endl << std::endl ;

    x_array[i] = x ;
    y_array[i] = y ;
    z_array[i] = z ;

    return  0 ;
}







PyObject *smc_parallel(PyObject *self, PyObject *args){
	PyObject *array = NULL ;
    PyObject *pList = NULL ;

    //std::ostringstream sstream;
    //std::string remark = "#hello";
    //const std::string filename = "dum.txt"; 

    //std::ofstream outfile(filename.c_str()) ;
    //outfile << remark << std::endl;
    int i, j, ki; 
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

    std::cout << "\n\nI am ready to return\n\n" << std::endl ;

    int accepted = 0 ;

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
       
        for(j = 0 ; j < natoms ; j++){ 
            accepted = surface_move(x_array, y_array, z_array, j) ;
        }

    //    std::cout << "c : x_array[0] = " << x_array[0] << std::endl ;
            // stepresult = dcdio.write_dcdstep(filepointer, tx, ty, tz, step)¬
        stepresult = write_dcdstep(filepointer, natoms, x_array, y_array, z_array, i) ;

    } // end of i loop

    stepresult = close_dcd_write(filepointer) ;

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

