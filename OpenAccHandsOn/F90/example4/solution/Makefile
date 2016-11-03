FC=pgfortran
FCFLAGS=-O3 -ta=tesla -Minfo=accel -Mcudalib=curand -Mcuda

all: main.F90
	${FC} ${FCFLAGS} main.F90 
clean:
	rm -Rf a.out pgprof*
