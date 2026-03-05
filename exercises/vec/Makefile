CXXFLAGS = -Wall -Wextra -std=c++23
LINK.0 = $(CXX) $(CXXFLAGS)
LDLIBS += -lstdc++ -lm 

Out.txt: main
	./main > Out.txt

main: main.o vec.o

main.o: main.cc vec.cc vec.h
# 	$(CXX) $(CXXFLAGS) -c main.cc

vec.o: vec.cc vec.h
# 	$(CXX) $(CXXFLAGS) -c vec.cc

.PHONY : clean
clean:
	$(RM) main main.o Out.txt *.o
