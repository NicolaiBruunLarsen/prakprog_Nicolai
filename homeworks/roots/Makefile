CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

OBJ = main.o vector.o matrix.o QR.o roots.o ODE.o

main: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o main

main.o: main.cc roots.h ODE.h vector.h matrix.h QR.h
	$(CXX) $(CXXFLAGS) -c main.cc

vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c vector.cc

matrix.o: matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

QR.o: QR.cc QR.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c QR.cc

roots.o: roots.cc roots.h QR.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c roots.cc

ODE.o: ODE.cc ODE.h vector.h
	$(CXX) $(CXXFLAGS) -c ODE.cc

Out.txt: main
	./main > Out.txt

run: Out.txt
	cat Out.txt

plot: Out.txt plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main *.o Out.txt out.*.txt fig_*.svg
