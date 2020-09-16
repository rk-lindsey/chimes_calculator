CC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 -fPIC

CHIMESFF_SRC=$(CC_LOC)/../../../../src/chimesFF.cpp
CHIMESFF_HDR=$(CC_LOC)/../../../../src/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)

WRAPPER_SRC=$(CC_LOC)/../../cpp_and_c/src/wrapper-C.cpp
WRAPPER_HDR=$(CC_LOC)/../../cpp_and_c/src/wrapper-C.h

wrapper-C.o : $(WRAPPER_SRC) $(WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(WRAPPER_SRC) $(CHIMESFF_SRC) -I $(CC_LOC)/../../cpp_and_c/src/ -I $(CC_LOC)/../../../../src/
	
wrapper-C.so: 	
	$(CXX) -shared -o libwrapper-C.so wrapper-C.o chimesFF.o


all:
	make chimesFF.o
	make wrapper-C.o
	make wrapper-C.so
	
clean:
	rm -f *.o *.so
