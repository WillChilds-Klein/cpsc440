CFLAGS = -g3 -std=c99 
CC = gcc

all: inv_double_gs

inv_double_gs: inv_double_gs.o
		${CC} ${CFLAGS} inv_double_gs.c -o inv_double_gs -lm