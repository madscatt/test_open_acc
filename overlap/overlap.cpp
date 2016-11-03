//#include <vector>
#include <iostream>
#include <math.h>  
#include <stdio.h>
#include <ctime>

/********* methods        ******************/

/********* main           ******************/

//void overlap(double x, double y, double z, int natoms) ;

void overlap(double x[], double y[], double z[], int natoms){

    float lcut= 2000.3 ;
    int check=0 ;
    for (int i=0; i< natoms-1 ; i++){
        float x1 = 35.8 ;
        float y1 = 25.8 ;
        float z1 = 75.8 ;
        for (int j=i+1; j<natoms ; j++){
            float x2 = 15.8 ;
            float y2 = 35.8 ;
            float z2 = 45.8 ;
            
            float sdist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) ;
            float dist=sqrt(sdist) ;
            if(dist<lcut){
                check=1 ;
                printf(" i = %d\tj = %d\tdist = %f\n", i,j,dist) ;
                printf(" x1 = %f\t y1 = %f\t z1 = %f\n", x1,y1,z1) ;
                printf(" x2 = %f\t y2 = %f\t z2 = %f\n", x2,y2,z2) ;
                break ;
                }
            }
            if(check==1) {
                break ;
            }
        }
    return ;
}


int main(){

	std::cout << "hello \n\n\n" ;
	
	//const std::string pdb_filename = "min3.pdb" ;

	//std::cout << "testing sasmol " << " " << pdb_filename << std::endl ; 

    //float rg = 38.4 ;
	//std::cout << "radius of gyration = " << rg << std::endl ;
	//std::cout << "sasmol.py : radius of gyration = 64.043168998442439 " << std::endl ;
	
	//clock_t tStart = clock();

	//float theta = acos(-1.0)/4.0 ;
	//std::cout << "theta = " << theta << std::endl ;

    //int natoms = 100 ;
    //double x[natoms]; double y[natoms]; double z[natoms] ;
    //for (int i = 0 ; i < natoms ; i++){
    //    x[i] = 0.0 ; y[i] = 0.0 ; z[i] = 0.0 ;
    //}
    //overlap(x, y, z, natoms) ;
	
	return 0 ;


}
