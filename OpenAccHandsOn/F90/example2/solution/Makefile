FC=pgfortran
FCFLAGS=-O3 -acc -ta=tesla -Minfo=accel

all: kernels.F90 main.F90
	${FC} ${FCFLAGS} kernels.F90 main.F90 
clean:
	rm -Rf a.out pgprof* *.o *.mod
