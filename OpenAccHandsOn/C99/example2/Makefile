CC=pgc++
CFLAGS=-fast -acc -Minfo -ta=tesla:cuda8

all:
	${CC} ${CFLAGS} main.cpp 
clean:
	rm -Rf a.out pgprof*
