CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

OBJECTS = vector.cc matrix.cc QR.cc min.cc ann.cc main.cc

main: $(OBJECTS) vector.h matrix.h QR.h min.h ann.h
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o main

Out.txt: main
	./main > Out.txt

run: Out.txt
	cat Out.txt

plot: Out.txt plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main Out.txt out.*.txt fig.*.svg *.o *.log
