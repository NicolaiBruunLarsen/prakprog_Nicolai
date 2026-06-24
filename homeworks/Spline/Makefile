CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

.PHONY: all run plot test clean

all: run plot

main: main.o Splines.o vector.o matrix.o
	$(CXX) $(CXXFLAGS) $^ -o $@

main.o: main.cc Splines.h vector.h
	$(CXX) $(CXXFLAGS) -c main.cc

Splines.o: Splines.cc Splines.h vector.h
	$(CXX) $(CXXFLAGS) -c Splines.cc

vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c vector.cc

matrix.o: matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

Out.txt: main
	./main > Out.txt

run: Out.txt
	cat Out.txt

plot: Out.txt plot.gpi
	gnuplot plot.gpi

test: main
	./main --test

clean:
	$(RM) main *.o data.txt data_lin.txt data_quad.txt Out.txt *.svg *.log