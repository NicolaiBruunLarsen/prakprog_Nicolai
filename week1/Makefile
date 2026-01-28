CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

OBJS = main.o hello.o
TARGET = main

all: $(TARGET)

Out.txt: $(TARGET)
	./$(TARGET) > Out.txt

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

hallo.o: main.cc hello.h
	$(CXX) $(CXXFLAGS) -c hallo.cc

clean:
	rm -f $(OBJS) $(TARGET)
