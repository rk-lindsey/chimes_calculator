
CXX   = g++ -O3 -std=c++11
DEBUG = 1

CXX_LOC = $(realpath .)
CXX_SRC = main.cpp $(CXX_LOC)/../../src/serial_chimes_interface.cpp $(CXX_LOC)/../../../chimesFF/src/chimesFF.cpp
CXX_HDR = $(CXX_LOC)/../../src/chimescalc_serial.h $(CXX_LOC)/../../../chimesFF/src/chimescalc.h

clean:

clean-all:
	rm -f CPP-interface

all: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o CPP-interface -I $(CXX_LOC)/../../src -I $(CXX_LOC)/../../../chimesFF/src/ -DDEBUG=${DEBUG}
