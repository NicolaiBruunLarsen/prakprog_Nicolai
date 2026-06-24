CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

all: run plot

main: main.cc vector.cc matrix.cc EVD.cc
	$(CXX) $(CXXFLAGS) main.cc vector.cc matrix.cc EVD.cc -o main

test: main
	./main --test

run: main
	./main -rmax 10 -dr 0.3

# Matrix sizes for scaling test
Ns := $(shell seq 20 10 140)

out.times.txt: main Makefile
	>$@
	for N in $(Ns); do \
		./main --time $$N >> $@; \
	done

# Parallel timing version
out.times.parallel.txt: main Makefile
	$(RM) log.*
	for N in $(Ns); do \
		./main --time $$N > log.$$N & \
	done; wait
	>$@
	for N in $(Ns); do cat log.$$N >> $@; done
	sort -n $@ -o $@
	$(RM) log.*

# Hydrogen calculation: energies and wavefunctions
out.wavefunctions.txt: main
	./main -rmax 10 -dr 0.3 > out.energies.txt

# Convergence with respect to dr
out.convergence.dr.txt: main
	./main --convergence-dr -rmax 10

# Convergence with respect to rmax
out.convergence.rmax.txt: main
	./main --convergence-rmax -dr 0.2

plot: out.times.txt out.times.parallel.txt out.wavefunctions.txt out.convergence.dr.txt out.convergence.rmax.txt plot_scaling.gpi plot.gpi
	gnuplot plot_scaling.gpi
	gnuplot plot.gpi

clean:
	$(RM) main \
	      out.times.txt out.times.parallel.txt \
	      out.wavefunctions.txt out.energies.txt \
	      out.convergence.dr.txt out.convergence.rmax.txt \
	      times.svg convergence.dr.svg convergence.rmax.svg wavefunctions.svg \
	      log.* fit.log *.o *.exe *.log
