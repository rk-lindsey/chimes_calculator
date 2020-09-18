
CXX=g++ -O3 -std=c++11

CXX_LOC = $(realpath .)
CXX_SRC = main.cpp $(CXX_LOC)/../../../serial_interface/src/serial_chimes_interface.cpp $(CXX_LOC)/../../src/chimesFF.cpp
CXX_HDR = $(CXX_LOC)/../../../serial_interface/src/serial_chimes_interface.h $(CXX_LOC)/../../src/chimesFF.h

test-CPP: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o test-CPP -I $(CXX_LOC)/../../src -I $(CXX_LOC)/../../../serial_interface/src/
clean:
	rm -f test-CPP
