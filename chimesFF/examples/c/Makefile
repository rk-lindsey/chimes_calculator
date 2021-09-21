CC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11

CHIMESFF_SRC=$(CC_LOC)/../../src/chimescalc.cpp
CHIMESFF_HDR=$(CC_LOC)/../../src/chimescalc.h

chimescalc.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)

WRAPPER_SRC=$(CC_LOC)/../../api/chimescalc_C.cpp
WRAPPER_HDR=$(CC_LOC)/../../api/chimescalc_C.h

chimescalc_C.o : $(WRAPPER_SRC) $(WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(WRAPPER_SRC) $(CHIMESFF_SRC) -I $(CC_LOC)/../../api/ -I $(CC_LOC)/../../src/

CC = gcc

test_wrapper-C.o : main.c $(WRAPPER_SRC) $(WRAPPER_HDR)
	$(CC) -c main.c -o test_wrapper-C.o -I $(CC_LOC)/../../api/

LINKS = chimescalc.o test_wrapper-C.o chimescalc_C.o

test_wrapper-C : $(LINKS)
	$(CXX)  $(LINKS) -o C_wrapper-direct_interface

clean:
	rm -f *.o

clean-all:
	rm -f *.o
	rm -f C_wrapper-direct_interface

all:
	make chimescalc.o
	make chimescalc_C.o
	make test_wrapper-C.o
	make test_wrapper-C
	make clean
