# Makefile for N-Body Simulation

CXX = g++
CXXFLAGS = -O2 -std=c++11
TARGET = nbody
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o