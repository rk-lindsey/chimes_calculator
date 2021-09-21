CC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 -fPIC

CHIMESFF_SRC=$(CC_LOC)/../../src/chimescalc.cpp
CHIMESFF_HDR=$(CC_LOC)/../../src/chimescalc.h

chimescalc.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)

CC_WRAPPER_SRC=$(CC_LOC)/../../api/chimescalc_C.cpp
CC_WRAPPER_HDR=$(CC_LOC)/../../api/chimescalc_C.h

chimescalc_C.o : $(CC_WRAPPER_SRC) $(CC_WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(CC_WRAPPER_SRC) $(CHIMESFF_SRC) -I $(CC_LOC)/../../api/ -I $(CC_LOC)/../../src/

chimescalc_C.so:
	$(CXX) -shared -o lib-C_wrapper-direct_interface.so chimescalc_C.o chimescalc.o


clean:
	rm -f *.o

clean-all:

	make clean
	rm -f *.so

all:
	make chimescalc.o
	make chimescalc_C.o
	make chimescalc_C.so
	make clean
