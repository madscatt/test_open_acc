#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;

//const unsigned int WIDTH=16384;
//const unsigned int HEIGHT=16384;
const unsigned int WIDTH=4096;
const unsigned int HEIGHT=4096;
const unsigned int MAX_ITERS=100;
const unsigned int MAX_COLOR=255;
const double xmin=-1.7;
const double xmax=.5;
const double ymin=-1.2;
const double ymax=1.2;
const double dx=(xmax-xmin)/WIDTH;
const double dy=(ymax-ymin)/HEIGHT;


unsigned char mandlebrot(int Px, int Py) {
  int i;
  double x0=xmin+Px*dx;
  double y0=ymin+Py*dy;
  double x=0.0;
  double y=0.0;
  for(i=0;i<MAX_ITERS;i++) {
    if (x*x+y*y>=4.0) break;
    double xtemp=x*x-y*y+x0;
    y=2*x*y+y0;
    x=xtemp;
  }
  return (double)MAX_COLOR*i/MAX_ITERS;
}

int main() {
  
  size_t bytes=WIDTH*HEIGHT*sizeof(unsigned char);
  unsigned char *image=(unsigned char*)malloc(bytes);
  FILE *fp=fopen("image.pgm","wb");
  fprintf(fp,"P5\n%s\n%d %d\n%d\n","#comment",WIDTH,HEIGHT,MAX_COLOR);

  #pragma acc parallel
  {
  #pragma acc loop gang
  for(int y=0;y<HEIGHT;y++) {
    #pragma acc loop vector
    for(int x=0;x<WIDTH;x++) {
      image[y*WIDTH+x]=mandlebrot(x,y);
    }
  }
  }
  
  fwrite(image,sizeof(unsigned char),WIDTH*HEIGHT,fp);
  fclose(fp);
  free(image);
  return 0;
}
