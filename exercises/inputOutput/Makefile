CXXFLAGS = -Wall -Wextra -std=c++23
LINK.0 = $(CXX) $(CXXFLAGS)
LDLIBS += -lstdc++ -lm 

Out.txt: main input.txt
	./main --input input.txt --output Out.txt

main: main.o 

main.o: main.cc input.txt
# 	$(CXX) $(CXXFLAGS) -c main.cc

.PHONY : clean
clean:
	$(RM) main main.o Out.txt *.o
