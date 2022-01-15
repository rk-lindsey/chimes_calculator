CXX     = mpicxx -O3 -std=c++11 -D USE_MPI
CCOMP   = mpicc -O3 -D USE_MPI
DEBUG   = 1

C_LOC = $(realpath .)


CHIMESFF_LOC=$(C_LOC)/../../../chimesFF/src
CHIMESFF_SRC=$(CHIMESFF_LOC)/chimesFF.cpp
CHIMESFF_HDR=$(CHIMESFF_LOC)/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC) -I $(CHIMESFF_LOC)


ser_LOC=$(C_LOC)/../../../serial_interface/src
ser_SRC=$(ser_LOC)/serial_chimes_interface.cpp
ser_HDR=$(mpi_LOC)/serial_chimes_interface.h


serial_chimes_interface.o : $(ser_SRC)
	$(CXX)  -c $(ser_SRC)  -I $(ser_LOC) -I $(CHIMESFF_LOC)

mpi_LOC=$(C_LOC)/../../src
mpi_SRC=$(mpi_LOC)/mpi_chimes_interface.cpp
mpi_HDR=$(mpi_LOC)/mpi_chimes_interface.h

mpi_chimes_interface.o : $(mpi_SRC)
	$(CXX) -c $(mpi_SRC) -c $(ser_SRC) -I $(mpi_LOC) -I $(ser_LOC) -I $(CHIMESFF_LOC)


mpi_WRAPPER_LOC=$(C_LOC)/../../api
mpi_WRAPPER_SRC=$(mpi_WRAPPER_LOC)/chimescalc_mpi_C.cpp
mpi_WRAPPER_HDR=$(mpi_WRAPPER_LOC)/chimescalc_mpi_C.h


chimescalc_mpi_C.o : $(mpi_WRAPPER_SRC)
	$(CXX) -c $(mpi_WRAPPER_SRC) -I $(mpi_WRAPPER_LOC) -I $(mpi_LOC)  -I $(CHIMESFF_LOC) -I $(ser_LOC) 

TEST_LOC=$(C_LOC)
TEST_SRC=$(TEST_LOC)/main.c


main.o: $(TEST_SRC)
	$(CCOMP) -c $(TEST_SRC) -DDEBUG=${DEBUG}  -I $(mpi_WRAPPER_LOC)   -I $(CHIMESFF_LOC) -I $(mpi_LOC) -I $(ser_LOC)  #   -I $(ser_WRAPPER_LOC) -I $(ser_LOC)

TEST_LNK = chimesFF.o serial_chimes_interface.o mpi_chimes_interface.o chimescalc_mpi_C.o main.o


test-C: $(TEST_LNK)
	$(CXX) $(TEST_LNK) -o chimescalc-test_mpi-C


clean:
	rm -f *.o


clean-all:
	make clean
	rm -f chimescalc-test_mpi-C


all:
	make chimesFF.o
	make serial_chimes_interface.o
	make mpi_chimes_interface.o
	make chimescalc_mpi_C.o
	make main.o
	make test-C
	make clean
