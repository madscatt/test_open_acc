CC=g++
FC=pgfortran
CFLAGS=-O3
FCFLAGS=-fast -Minfo=all

all: main.cpp matrix.F90 vector.F90 main.F90
	${CC} ${CFLAGS} main.cpp -o cmain
	${FC} ${FCFLAGS} matrix.F90 vector.F90 main.F90 -o fmain
clean:
	rm -Rf a.out pgprof* cmain fmain *.o *.mod
