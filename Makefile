CFLAGS = -g3 -std=c99 
CC = gcc

all: upperhes

upperhes: upperhes.o
		${CC} ${CFLAGS} upperhes.c -o upperhes -lm