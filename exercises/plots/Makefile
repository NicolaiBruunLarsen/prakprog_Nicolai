CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra
LDLIBS = -lm

all: plot.erf.svg plot.gamma.svg plot.lngamma.svg

main: main.o sfuns.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

main.o: main.cc sfuns.h
	$(CXX) $(CXXFLAGS) -c main.cc

sfuns.o: sfuns.cc sfuns.h
	$(CXX) $(CXXFLAGS) -c sfuns.cc

data.txt: main
	./main > $@

plot.erf.svg: data.txt values.dat Makefile
	echo '\
	set terminal svg background "white";\
	set output "$@";\
	set key left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" index 0 with lines title "erf(x)"\
	,"values.dat" using 1:2 with points pointtype 4 title "tabulated erf(x)"\
	' | gnuplot

plot.gamma.svg: data.txt Makefile
	echo '\
	set terminal svg background "white";\
	set output "$@";\
	set key bottom;\
	set xlabel "x";\
	set ylabel "Γ(x)";\
	set title "Gamma function";\
	plot [-4:7][-20:130] \
	"$<" index 1 with lines title "gamma(x)"\
	,"$<" index 3 with points pointtype 4 title "factorials"\
	' | gnuplot

plot.lngamma.svg: data.txt Makefile
	echo '\
	set terminal svg background "white";\
	set output "$@";\
	set key left;\
	set xlabel "x";\
	set ylabel "ln(Γ(x))";\
	set title "Log Gamma function";\
	plot [0:12][-3:12]\
	"$<" index 2 with lines title "ln(gamma(x))"\
	' | gnuplot

clean:
	$(RM) main *.o data.txt *.svg