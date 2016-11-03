FC=pgfortran
FCFLAGS=-fast -Minfo=accel -O3 -acc -ta=nvidia

a.out: matrix.o vector.o main.F90
	${FC} ${FCFLAGS} main.F90 vector.o matrix.o

matrix.o: matrix.F90
	${FC} ${FCFLAGS} matrix.F90 -c

vector.o: vector.F90
	${FC} ${FCFLAGS} vector.F90 -c

clean:
	rm -Rf a.out pgprof* cmain fmain *.o *.mod
