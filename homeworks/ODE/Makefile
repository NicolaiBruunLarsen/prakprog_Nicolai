CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

main: main.cc ODE.cc vector.cc
	$(CXX) $(CXXFLAGS) main.cc ODE.cc vector.cc -o main

run: main
	./main

plot: run plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main out.*.txt fig_*.svg *.o *.exe *.log
