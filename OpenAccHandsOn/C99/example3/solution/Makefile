CC=pgc++
CFLAGS=-acc -Minfo=accel -ta=nvidia -O3

all:
	${CC} ${CFLAGS} main.cpp 
clean:
	rm -Rf a.out pgprof* image.pgm
