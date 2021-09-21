
CXX=g++ -O3 -std=c++11

CXX_LOC = $(realpath .)
CXX_SRC = main.cpp $(CXX_LOC)/../../../serial_interface/src/chimescalc_serial.cpp $(CXX_LOC)/../../src/chimescalc.cpp
CXX_HDR = $(CXX_LOC)/../../../serial_interface/src/chimescalc_serial.h $(CXX_LOC)/../../src/chimescalc.h

all: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o CPP-interface -I $(CXX_LOC)/../../src -I $(CXX_LOC)/../../../serial_interface/src/
clean:
	rm -f CPP-interface
