CFLAGS = -g3 -std=c99 
CC = gcc

all: qr_symmetric

qr_symmetric: qr_symmetric.o
		${CC} ${CFLAGS} qr_symmetric.c -o qr_symmetric -lm