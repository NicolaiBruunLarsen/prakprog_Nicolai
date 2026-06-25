CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

all: run plot



main: main.cc ann.cc ann.h min.cc min.h QR.cc QR.h matrix.cc matrix.h vector.cc vector.h
	$(CXX) $(CXXFLAGS) main.cc ann.cc min.cc QR.cc matrix.cc vector.cc -o main

Out.txt: main
	./main > Out.txt

run: Out.txt
	cat Out.txt

plot: Out.txt plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main Out.txt out.*.txt fig.*.svg *.o *.log
