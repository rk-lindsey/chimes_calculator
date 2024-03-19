FCC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 -fPIC

CHIMESFF_SRC=$(FCC_LOC)/../../src/chimesFF.cpp
CHIMESFF_HDR=$(FCC_LOC)/../../src/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)
	
CC_WRAPPER_SRC=$(FCC_LOC)/../../api/chimescalc_C.cpp
CC_WRAPPER_HDR=$(FCC_LOC)/../../api/chimescalc_C.h
		
chimescalc_C.o : $(CC_WRAPPER_SRC) $(CC_WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(CC_WRAPPER_SRC) $(CHIMESFF_SRC) -I $(FCC_LOC)/../../api/ -I $(FCC_LOC)/../../src/	

FCC = gfortran -O3 -fPIC -std=f2003

FCC_WRAPPER_SRC=$(FCC_LOC)/../../api/chimescalc_F.f90
	
chimescalc_F.o chimescalc.mod : $(FCC_WRAPPER_SRC)
	$(FCC) -c $(FCC_WRAPPER_SRC) -o chimescalc_F.o
	
main.o : main.F90
	$(FCC) -c main.F90 -o main.o

LINKS = chimesFF.o chimescalc_F.o main.o chimescalc_C.o

# Rudimentary OS detection
UNAME := $(shell uname)

ifeq (${UNAME},Darwin)
	FCC += -lc++
else
	FCC += -lstdc++
endif

test_wrapper-F: $(LINKS)
	$(FCC) $(LINKS) -o chimescalc-test_direct-F

clean:
	rm -f *.o *.mod

clean-all:
	rm -f *.o *.mod
	rm -f chimescalc-test_direct-F

all:
	make chimesFF.o
	make chimescalc_C.o
	make chimescalc_F.o chimescalc.mod
	make main.o 
	make test_wrapper-F
	#make clean
