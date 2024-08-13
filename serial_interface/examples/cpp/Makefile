
CXX   = nvc++ -O3 -std=c++17
DEBUG = 1

CXX_LOC = $(realpath .)
CXX_SRC = main.cpp $(CXX_LOC)/../../src/serial_chimes_interface.cpp $(CXX_LOC)/../../../chimesFF/src/chimesFF.cpp
CXX_HDR = $(CXX_LOC)/../../src/serial_chimes_interface.h $(CXX_LOC)/../../../chimesFF/src/chimesFF.h

clean:

clean-all:
	rm -f chimescalc

all: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o chimescalc -I $(CXX_LOC)/../../src -I $(CXX_LOC)/../../../chimesFF/src/ -DDEBUG=${DEBUG}
