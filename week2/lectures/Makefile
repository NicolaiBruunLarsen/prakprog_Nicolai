CXXFLAGS = -Wall -Wextra -std=c++23
LINK.0 = $(CXX) $(CXXFLAGS)
LDLIBS += -lstdc++ -lm 

Out .txt : main
	./main > Out.txt

main : main.o hello.o sfuns.o
# 	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o main
# This line already build in make

main.o : hello.h sfuns.h main.cc

hello.o : hello.h hello.cc

sfuns.o : sfuns.h sfuns.cc

.PHONY : clean
clean :
	$(RM) main main.o Out.txt *.o








# CXX = g++
# CXXFLAGS = -Wall -Wextra -std=c++17

# OBJS = main.o hello.o
# TARGET = main

# all: $(TARGET)

# Out.txt: $(TARGET)
# 	./$(TARGET) > Out.txt

# $(TARGET): $(OBJS)
# 	$(CXX) $(OBJS) -o $(TARGET)

# hallo.o: main.cc hello.h
# 	$(CXX) $(CXXFLAGS) -c hallo.cc

# clean:
# 	rm -f $(OBJS) $(TARGET)
