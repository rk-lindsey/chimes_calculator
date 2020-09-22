
CXX   = g++ -O3 -std=c++11
DEBUG = 0

CXX_LOC = $(realpath .)
CXX_SRC = main.cpp $(CXX_LOC)/../../src/serial_chimes_interface.cpp $(CXX_LOC)/../../../chimesFF/src/chimesFF.cpp
CXX_HDR = $(CXX_LOC)/../../src/serial_chimes_interface.h $(CXX_LOC)/../../../chimesFF/src/chimesFF.h


clean-all:
	rm -f test-CPP

all: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o test-CPP -I $(CXX_LOC)/../../src -I $(CXX_LOC)/../../../chimesFF/src/ -DDEBUG=${DEBUG}
