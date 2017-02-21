CC=/share/apps/local/pgi/linux86-64/2016/bin/pgc++
#CC=/opt/pgi/osx86-64/16.10/bin/pgc++
#CC=/opt/pgi/osx86/16.10/bin/pgc++
#CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,fastmath -O3
#CFLAGS= -c -fast -acc -Minfo=accel -ta=tesla:cuda8,maxregcount:63 -O3
## to keep cuda code
## open_acc_tests_gv.n001.gpu <- has the CUDA code ... not compatiable with 
#CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,nollvm,keepgpu -O3 -v 
#IFLAGS= /opt/pgi/osx86-64/16.10/include
#CFLAGS= -c -g -Mlist -Minfo=accel -ta=multicore -O3
#CFLAGS= -fPIC -c -Minfo=accel -ta=multicore -O3 -shared 
#CFLAGS= -fPIC -c -Minfo=all -ta=multicore -O3 -shared 
CFLAGS= -fPIC -c -fast -acc -Minfo=accel -ta=tesla:cuda8,fastmath -O3 -shared
CFLAGS= -fPIC -c -fast -acc -Minfo=accel -ta=tesla:nordc,fastmath -O3 -shared
#CFLAGS= -fPIC -c -fast -acc -Minfo=accel -ta=multicore -O3 -shared
#CFLAGS= -fPIC -c -fast -acc -Minfo=accel -ta=multicore -O3 -shared 

#CC=g++
#CFLAGS= -fPIC -c -O3 -shared 

#gcc -c -Wall -Werror -fPIC shared.c
#gcc -shared -o libshared.so shared.o
# 
# NOTE: shared object file (linux) needs to be /share/apps/local/anaconda/lib
# 
DEST=/share/apps/local/anaconda/lib

all:
	${CC} ${CFLAGS} oacc_pr.cpp 
	${CC} -shared -o liboacc_pr.so oacc_pr.o
	ar -cvq liboacc_pr.a oacc_pr.o
clean:
	rm -Rf *.o *.so pgprof*
install:
	cp liboacc_pr.so liboacc_pr.a ${DEST}
