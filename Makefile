CC = g++
CFLAGS = -Wall
DEBUG = -g 
OPT = -O3 -fopenmp
LIBS = ~/gsl/libs/libgsl.a
INC = -I/usr/lib

OBJS = main.o chain.o graph.o
TARGET = IRv1

$(TARGET): $(OBJS)
	$(CC)	$(OPT)	$(INC)	-o	$(TARGET) -fopenmp -m64  main.cpp chain.cpp graph.cpp	$(LIBS)

clean: 
	-rm	-f	$(TARGET)	$(OBJS)