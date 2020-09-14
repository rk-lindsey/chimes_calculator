
CXX=g++ -O3 -std=c++11

CXX_LOC = $(realpath .)
CXX_SRC = main.cpp ../src/serial_chimes_interface.cpp ../../../src/chimesFF.cpp
CXX_HDR = $(CXX_LOC)/../src/serial_chimes_interface.h $(CXX_LOC)/../../../src/chimesFF.h

example: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o example -I ../src -I ../../../src
clean:
	rm -f example
