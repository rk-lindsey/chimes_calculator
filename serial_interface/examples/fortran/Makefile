CXX     = g++ -O3 -std=c++11 -fPIC
FCC     = gfortran -O3 -fPIC -std=f2003 -cpp # last flag allows processing of C-style pre-processor directives

DEBUG   = 1

#CCFLAGS = -fPIC -O3 -std=c++11
#FFLAGS  = -O3 -fPIC
#COMM1   = chimescalc_serial_F.o main.o

F_LOC = $(realpath .)

chimescalc_LOC=$(F_LOC)/../../../chimesFF/src
chimescalc_SRC=$(chimescalc_LOC)/chimescalc.cpp
chimescalc_HDR=$(chimescalc_LOC)/chimescalc.h

chimescalc.o : $(chimescalc_SRC)
	$(CXX) -c $(chimescalc_SRC) -I $(chimescalc_LOC)

SERIAL_LOC=$(F_LOC)/../../src
SERIAL_SRC=$(SERIAL_LOC)/chimescalc_serial.cpp
SERIAL_HDR=$(SERIAL_LOC)/chimescalc_serial.h

chimescalc_serial.o : $(SERIAL_SRC)
	$(CXX) -c $(SERIAL_SRC) -I $(SERIAL_LOC) -I $(chimescalc_LOC)

WRAPPERC_LOC=$(F_LOC)/../../api
WRAPPERC_SRC=$(WRAPPERC_LOC)/chimescalc_serial_C.cpp
WRAPPERC_HDR=$(WRAPPERC_LOC)/wrapper-h.cpp

chimescalc_serial_C.o : $(WRAPPERC_SRC)
	$(CXX) -c $(WRAPPERC_SRC) -I $(WRAPPERC_LOC) -I $(SERIAL_LOC) -I $(chimescalc_LOC)

WRAPPERF_LOC=$(F_LOC)/../../api
WRAPPERF_SRC=$(WRAPPERF_LOC)/chimescalc_serial_F.F90

chimescalc_serial_F.o chimescalc_serial.mod : $(WRAPPERF_SRC)
	$(FCC) -c $(WRAPPERF_SRC) -o chimescalc_serial_F.o

main.o : main.F90
	$(FCC) -c main.F90 -o main.o  -DDEBUG=${DEBUG}


LINKS = chimescalc.o chimescalc_serial.o  chimescalc_serial_C.o chimescalc_serial_F.o main.o

# Rudimentary OS detection
UNAME := $(shell uname)

ifeq (${UNAME},Darwin)
	FCC += -lc++
else
	FCC += -lstdc++
endif

test-F: $(LINKS)
	$(FCC) $(LINKS) -o fortran_wrapper-serial_interface
lib : chimescalc.o chimescalc_serial.o chimescalc_serial_C.o
	ar rcs libchimes.a chimescalc.o chimescalc_serial.o chimescalc_serial_C.o

clean:
	rm -f *.o *.mod

clean-all:
	make clean
	rm -f fortran_wrapper-serial_interface
	rm -f libchimes.a

all:
	make chimescalc.o
	make chimescalc_serial.o
	make chimescalc_serial_C.o
	make chimescalc_serial_F.o chimescalc_serial.mod
	make main.o
	make test-F
	make lib
	make clean
