CFLAGS = -g3 -std=c99 -w
CC = gcc

all: inv_double_gs upperhes qr_symmetric

inv_double_gs: inv_double_gs.o
	${CC} ${CFLAGS} inv_double_gs.c -o inv_double_gs -lm

upperhes: upperhes.o
	${CC} ${CFLAGS} upperhes.c -o upperhes -lm

qr_symmetric: qr_symmetric.o
	${CC} ${CFLAGS} qr_symmetric.c -o qr_symmetric -lm

clean:
	rm -rf *.o *.dSYM inv_double_gs upperhes qr_symmetric