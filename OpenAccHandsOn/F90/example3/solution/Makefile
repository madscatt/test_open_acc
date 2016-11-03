FC=pgfortran
FCFLAGS=-acc -Minfo=accel -ta=nvidia -O3

all: main.F90
	${FC} ${FCFLAGS} main.F90 
clean:
	rm -Rf a.out pgprof* image.pgm *.o *.mod
