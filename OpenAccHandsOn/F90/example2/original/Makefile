FC=gfortran
FCFLAGS=-g

all: kernels.F90 main.F90
	${FC} ${FCFLAGS} kernels.F90 main.F90 
clean:
	rm -Rf a.out pgprof* *.o *.mod
