CXX     = g++ -O3 -std=c++11
CCOMP   = gcc -O3
DEBUG   = 1

C_LOC = $(realpath .)

CHIMESFF_LOC=$(C_LOC)/../../../chimesFF/src
CHIMESFF_SRC=$(CHIMESFF_LOC)/chimescalc.cpp
CHIMESFF_HDR=$(CHIMESFF_LOC)/chimescalc.h

chimescalc.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC) -I $(CHIMESFF_LOC)

SERIAL_LOC=$(C_LOC)/../../src
SERIAL_SRC=$(SERIAL_LOC)/chimescalc_serial.cpp
SERIAL_HDR=$(SERIAL_LOC)/chimescalc_serial.h

chimescalc_serial.o : $(SERIAL_SRC)
	$(CXX) -c $(SERIAL_SRC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC)

WRAPPER_LOC=$(C_LOC)/../../api
WRAPPER_SRC=$(WRAPPER_LOC)/chimescalc_serial_C.cpp
WRAPPER_HDR=$(WRAPPER_LOC)/chimescalc_serial_C.h

chimescalc_serial_C.o : $(WRAPPER_SRC)
	$(CXX) -c $(WRAPPER_SRC) -I $(WRAPPER_LOC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC)

TEST_LOC=$(C_LOC)
TEST_SRC=$(TEST_LOC)/main.c

main.o: $(TEST_SRC)
	$(CCOMP) -c $(TEST_SRC) -DDEBUG=${DEBUG}  -I $(WRAPPER_LOC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC)

TEST_LNK = chimescalc.o chimescalc_serial.o chimescalc_serial_C.o main.o

test-C: $(TEST_LNK)
	$(CXX) $(TEST_LNK) -o C_wrapper-serial_interface

clean:
	rm -f *.o

clean-all:
	make clean
	rm -f C_wrapper-serial_interface

all:
	make chimescalc.o
	make chimescalc_serial.o
	make chimescalc_serial_C.o
	make main.o
	make test-C
	make clean
