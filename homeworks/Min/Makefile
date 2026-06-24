CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

OBJ = main.o vector.o matrix.o QR.o min.o

all: run plot

main: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o main

main.o: main.cc min.h QR.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c main.cc

vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c vector.cc

matrix.o: matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

QR.o: QR.cc QR.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c QR.cc

min.o: min.cc min.h QR.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c min.cc

Out.txt: main higgs.data.txt
	./main < higgs.data.txt > Out.txt 2> log

run: Out.txt
	cat Out.txt

plot: Out.txt plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main *.o Out.txt log out.*.txt fig.*.svg fit.log
