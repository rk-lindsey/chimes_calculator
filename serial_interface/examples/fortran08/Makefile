CXX     = g++ -O3 -std=c++11 -fPIC
FCC     = gfortran -O3 -fPIC -std=f2008 -cpp # last flag allows processing of C-style pre-processor directives

DEBUG   = 1

#CCFLAGS = -fPIC -O3 -std=c++11
#FFLAGS  = -O3 -fPIC
#COMM1   = chimescalc_serial_F.o main.o

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
WRAPPERC_SRC=$(WRAPPERC_LOC)/chimescalc_serial_C.cpp
WRAPPERC_HDR=$(WRAPPERC_LOC)/chimescalc_serial_.h 

chimescalc_serial_C.o : $(WRAPPERC_SRC)
	$(CXX) -c $(WRAPPERC_SRC) -I $(WRAPPERC_LOC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC)

WRAPPERF_LOC=$(F_LOC)/../../api
WRAPPERF_SRC=$(WRAPPERF_LOC)/chimescalc_serial_F.f90
WRAPPERF08_SRC=$(WRAPPERF_LOC)/chimescalc_serial_F08.f90

chimescalc_serial_F.o wrapper.mod : $(WRAPPERF_SRC)
	$(FCC) -c $(WRAPPERF_SRC) -o chimescalc_serial_F.o

chimescalc_serial_F08.o chimes_serial08.mod : $(WRAPPERF08_SRC)
	$(FCC) -c $(WRAPPERF08_SRC) -o chimescalc_serial_F08.o

main.o : main.F90
	$(FCC) -c main.F90 -o main.o  -DDEBUG=${DEBUG}


LINKS = chimesFF.o serial_chimes_interface.o  chimescalc_serial_C.o chimescalc_serial_F.o chimescalc_serial_F08.o main.o

# Rudimentary OS detection
UNAME := $(shell uname)

ifeq (${UNAME},Darwin)
	FCC += -lc++
else
	FCC += -stdc++
endif

test-F: $(LINKS)
	$(FCC) $(LINKS) -o fortran08_wrapper-serial_interface
lib : chimesFF.o serial_chimes_interface.o chimescalc_serial_C.o
	ar rcs libchimes.a chimesFF.o serial_chimes_interface.o chimescalc_serial_C.o

clean:
	rm -f *.o *.mod

clean-all:
	make clean
	rm -f fortran08_wrapper-serial_interface
	rm -f libchimes.a

all:
	make chimesFF.o
	make serial_chimes_interface.o
	make chimescalc_serial_C.o
	make chimescalc_serial_F.o chimescalc_serial.mod
	make chimescalc_serial_F08.o chimes_serial08.mod
	make main.o
	make test-F
	make lib
	make clean
