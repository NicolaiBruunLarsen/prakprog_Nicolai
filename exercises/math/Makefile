CXXFLAGS = -Wall -Wextra -std=c++23
LINK.0 = $(CXX) $(CXXFLAGS)
LDLIBS += -lstdc++ -lm 

Out.txt: main
	./main > Out.txt

main: main.o sfuns.o
# 	$(CXX) $(CXXFLAGS) -o  main

main.o: main.cc sfuns.h
# 	$(CXX) $(CXXFLAGS) -c main.cc

sfuns.o: sfuns.cc sfuns.h
# 	$(CXX) $(CXXFLAGS) -c sfuns.cc

.PHONY : clean
clean:
	$(RM) main main.o Out.txt *.o
