FCC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 -fPIC

CHIMESFF_SRC=$(FCC_LOC)/../../src/chimesFF.cpp
CHIMESFF_HDR=$(FCC_LOC)/../../src/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)
	
CC_WRAPPER_SRC=$(FCC_LOC)/../../api/wrapper-C.cpp
CC_WRAPPER_HDR=$(FCC_LOC)/../../api/wrapper-C.h
		
wrapper-C.o : $(CC_WRAPPER_SRC) $(CC_WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(CC_WRAPPER_SRC) $(CHIMESFF_SRC) -I $(FCC_LOC)/../../api/ -I $(FCC_LOC)/../../src/	

FCC = gfortran -O3 -fPIC -std=f2003

FCC_WRAPPER_SRC=$(FCC_LOC)/../../api/wrapper-F.F90
	
wrapper-F.o chimes.mod : $(FCC_WRAPPER_SRC)
	$(FCC) -c $(FCC_WRAPPER_SRC) -o wrapper-F.o	
	
main.o : main.F90
	$(FCC) -c main.F90 -o main.o

LINKS = chimesFF.o wrapper-F.o main.o wrapper-C.o

# Rudimentary OS detection
UNAME := $(shell uname)

ifeq (${UNAME},Darwin)
	FCC += -lc++
else
	FCC += -stdc++
endif

test_wrapper-F: $(LINKS)
	$(FCC) $(LINKS) -o test_wrapper-F

clean:
	rm -f *.o *.mod

clean-all:
	rm -f *.o *.mod
	rm -f test_wrapper-F

all:
	make chimesFF.o
	make wrapper-C.o
	make wrapper-F.o chimes.mod
	make main.o 
	make test_wrapper-F
	make clean
