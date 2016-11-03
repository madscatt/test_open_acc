CC=/share/apps/local/pgi/linux86-64/2016/bin/pgc++
CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8 -O3
#CFLAGS=-Minfo=accel -ta=host -O3

all:
	${CC} ${CFLAGS} main.cpp
clean:
	rm -Rf a.out pgprof*
