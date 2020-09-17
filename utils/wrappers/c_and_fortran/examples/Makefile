FCC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 -fPIC

CHIMESFF_SRC=$(FCC_LOC)/../../../../src/chimesFF.cpp
CHIMESFF_HDR=$(FCC_LOC)/../../../../src/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)
	
CC_WRAPPER_SRC=$(FCC_LOC)/../../cpp_and_c/src/wrapper-C.cpp
CC_WRAPPER_HDR=$(FCC_LOC)/../../cpp_and_c/src/wrapper-C.h
		
wrapper-C.o : $(CC_WRAPPER_SRC) $(CC_WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(CC_WRAPPER_SRC) $(CHIMESFF_SRC) -I $(FCC_LOC)/../../cpp_and_c/src/ -I $(FCC_LOC)/../../../../src/	

FCC = gfortran -O3 -fPIC

FCC_WRAPPER_SRC=$(FCC_LOC)/../src/wrapper-F.F90
	
wrapper-F.o chimes.mod : $(FCC_WRAPPER_SRC)
	$(FCC) -c $(FCC_WRAPPER_SRC) -o wrapper-F.o	
	
test_chimes_c_and_fortran_wrapper.o : test_chimes_c_and_fortran_wrapper.F90
	$(FCC) -c test_chimes_c_and_fortran_wrapper.F90 -o test_chimes_c_and_fortran_wrapper.o

LINKS = chimesFF.o wrapper-F.o test_chimes_F_api.o wrapper-C.o 

# Rudimentary OS detection
UNAME := $(shell uname)

ifeq (${UNAME},Darwin)
	FCC += -lc++
else
	FCC += -stdc++
endif

test_chimes_c_and_fortran_wrapper: $(LINKS)
	$(FCC) $(LINKS) -o test_chimes_c_and_fortran_wrapper.x

all:
	make chimesFF.o
	make wrapper-C.o
	make wrapper-F.o chimes.mod
	make test_chimes_c_and_fortran_wrapper.o 
	test_chimes_c_and_fortran_wrapper

clean:
	rm -f *.o *.mod

clean-all:
	rm -f *.o *.mod
	rm -f test_chimes_c_and_fortran_wrapper
