0 5.450000
1 32.020000
2 16.790000
3 2.910000 3.330000

CC=/share/apps/local/pgi/linux86-64/2016/bin/pgc++
CFLAGS= -fast -acc -Minfo=accel -ta=tesla:cuda8 -O3

nvprof -o timeline_3.out ./a.out
