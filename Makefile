all: Lab2

CXX = g++ -std=c++11

LINKFLAGS = -pedantic -Wall -fomit-frame-pointer -funroll-all-loops -O3 -lm

Lab2:
	$(CXX) $(LINKFLAGS) 312510155.cpp -o 312510155

clean:
	rm -rf *.o *.gch




