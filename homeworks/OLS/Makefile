CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra
LDLIBS = -lm

main: main.cc vector.cc matrix.cc QR.cc lsfit.cc
	$(CXX) $(CXXFLAGS) main.cc vector.cc matrix.cc QR.cc lsfit.cc -o main $(LDLIBS)

test: main
	./main --test

run: main
	./main --fit > out.results.txt

out.data.txt out.fit.txt out.logfit.txt out.results.txt: main
	./main --fit > out.results.txt

plot: out.data.txt out.fit.txt out.logfit.txt plot.gpi
	gnuplot plot.gpi

clean:
	$(RM) main out.data.txt out.fit.txt out.logfit.txt out.results.txt plot.fit.svg plot.logfit.svg log.* *.o *.exe *.log
