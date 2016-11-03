CC=pgc++
CFLAGS=-acc -ta=tesla -O3 -Minfo=accel

all: 
	${CC} ${CFLAGS} main.cpp 
clean:
	rm -Rf a.out pgprof*
