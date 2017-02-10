CC=/share/apps/local/pgi/linux86-64/2016/bin/pgc++
#CC=/opt/pgi/osx86-64/16.10/bin/pgc++
#CC=/opt/pgi/osx86/16.10/bin/pgc++
CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,fastmath -O3
CFLAGS= -c -fast -acc -Minfo=accel -ta=tesla:cuda8,maxregcount:63 -O3
## to keep cuda code
## open_acc_tests_gv.n001.gpu <- has the CUDA code ... not compatiable with 
#CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,nollvm,keepgpu -O3 -v 
#IFLAGS= /opt/pgi/osx86-64/16.10/include
#CFLAGS= -c -g -Mlist -Minfo=accel -ta=multicore -O3
CFLAGS= -fPIC -c -Minfo=accel -ta=multicore -O3 -shared 

all:
	${CC} ${CFLAGS} oacc_pr.cpp 
clean:
	rm -Rf *.so pgprof*
