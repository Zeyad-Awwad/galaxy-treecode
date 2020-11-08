CC=g++
CPPFLAGS=-std=c++11 -Wall -pedantic -O3
LDFLAGS= -O3

all: main

main: main.o Point.o Node.o Octree.o

Point.o: Point.hpp

Node.o: Point.hpp Node.hpp

Octree.o: Point.hpp Node.hpp Octree.hpp

main.o: main.cpp Point.hpp Node.hpp Octree.hpp

clean:
	rm -f all *.o core*