CXX     = g++ -O3 -std=c++11
CCOMP   = gcc -O3
DEBUG   = 1

C_LOC = $(realpath .)

CHIMESFF_LOC=$(C_LOC)/../../../chimesFF/src
CHIMESFF_SRC=$(CHIMESFF_LOC)/chimesFF.cpp
CHIMESFF_HDR=$(CHIMESFF_LOC)/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC) -I $(CHIMESFF_LOC)

SERIAL_LOC=$(C_LOC)/../../src
SERIAL_SRC=$(SERIAL_LOC)/serial_chimes_interface.cpp
SERIAL_HDR=$(SERIAL_LOC)/serial_chimes_interface.h

serial_chimes_interface.o : $(SERIAL_SRC)
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

TEST_LNK = chimesFF.o serial_chimes_interface.o chimescalc_serial_C.o main.o

test-C: $(TEST_LNK)
	$(CXX) $(TEST_LNK) -o chimescalc-test_serial-C

clean:
	rm -f *.o

clean-all:
	make clean
	rm -f chimescalc-test_serial-C

all:
	make chimesFF.o
	make serial_chimes_interface.o
	make chimescalc_serial_C.o
	make main.o
	make test-C
	make clean
