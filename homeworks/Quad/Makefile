CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

all: run plot

main: main.cc QUAD.cc QUAD.h
	$(CXX) $(CXXFLAGS) main.cc QUAD.cc -o main

Out.txt: main
	./main > Out.txt

run: Out.txt
	cat Out.txt

plot: Out.txt plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main Out.txt out.*.txt fig_*.svg *.o *.log fit.log
