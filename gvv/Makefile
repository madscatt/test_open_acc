CC=/share/apps/local/pgi/linux86-64/2016/bin/pgc++
CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,fastmath -O3
CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,maxregcount:63 -O3
## to keep cuda code
## open_acc_tests_gv.n001.gpu <- has the CUDA code ... not compatiable with 
#CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8,nollvm,keepgpu -O3 -v 
#CFLAGS=-Minfo=accel -ta=host -O3

all:
	${CC} ${CFLAGS} open_acc_tests_gv.cpp
clean:
	rm -Rf a.out pgprof*
