CXX     = g++ -g -std=c++11 -fPIC
FCC     = gfortran -g -fPIC -std=f2003 -cpp # last flag allows processing of C-style pre-processor directives

DEBUG   = 0

#CCFLAGS = -fPIC -O3 -std=c++11 
#FFLAGS  = -O3 -fPIC
#COMM1   = wrapper-F.o main.o 




F_LOC = $(realpath .)

CHIMESFF_LOC=$(F_LOC)/../../../chimesFF/src
CHIMESFF_SRC=$(CHIMESFF_LOC)/chimesFF.cpp
CHIMESFF_HDR=$(CHIMESFF_LOC)/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC) -I $(CHIMESFF_LOC)

SERIAL_LOC=$(F_LOC)/../../src
SERIAL_SRC=$(SERIAL_LOC)/serial_chimes_interface.cpp
SERIAL_HDR=$(SERIAL_LOC)/serial_chimes_interface.h

serial_chimes_interface.o : $(SERIAL_SRC) 
	$(CXX) -c $(SERIAL_SRC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC) 
	
WRAPPERC_LOC=$(F_LOC)/../../api
WRAPPERC_SRC=$(WRAPPERC_LOC)/wrapper-C.cpp
WRAPPERC_HDR=$(WRAPPERC_LOC)/wrapper-h.cpp

wrapper-C.o : $(WRAPPERC_SRC) 
	$(CXX) -c $(WRAPPERC_SRC) -I $(WRAPPERC_LOC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC)	

WRAPPERF_LOC=$(F_LOC)/../../api
WRAPPERF_SRC=$(WRAPPERF_LOC)/wrapper-F.F90

wrapper-F.o chimes.mod : $(WRAPPERF_SRC) 
	$(FCC) -c $(WRAPPERF_SRC) -o wrapper-F.o

main.o : main.F90
	$(FCC) -c main.F90 -o main.o  -DDEBUG=${DEBUG}


LINKS = chimesFF.o serial_chimes_interface.o  wrapper-C.o wrapper-F.o main.o

# Rudimentary OS detection
UNAME := $(shell uname)

ifeq (${UNAME},Darwin)
	FCC += -lc++
else
	FCC += -stdc++
endif

test-F: $(LINKS)
	$(FCC) $(LINKS) -o test-F
lib : chimesFF.o serial_chimes_interface.o wrapper-C.o
	ar rcs libchimes.a chimesFF.o serial_chimes_interface.o wrapper-C.o 

clean: 
	rm -f *.o *.mod
	
clean-all:
	make clean
	rm -f test-F
	
all:
	make chimesFF.o
	make serial_chimes_interface.o 
	make wrapper-C.o
	make wrapper-F.o chimes.mod
	make main.o
	make test-F
	make lib



