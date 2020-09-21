FCC = gfortran
CC=g++
CCFLAGS = -fPIC -O3 -std=c++11 
FFLAGS = -O3 -fPIC
COMM1 = wrapper-F.o main.o 

all:	lib test_F 

test_F : $(COMM1) 
	$(FCC) $(FFLAGS) $(COMM1) -o test_F.x libchimes.a -lc++ 

chimesFF.o : chimesFF.cpp
	$(CC) $(CCFLAGS) -c chimesFF.cpp

serial_chimes_interface.o : serial_chimes_interface.cpp
	$(CC) $(CCFLAGS) -c serial_chimes_interface.cpp

main.o : main.F90
	$(FCC) $(FFLAGS) -c main.F90

wrapper-C.o : wrapper-C.cpp
	$(CC) $(CCFLAGS) -c wrapper-C.cpp

lib : chimesFF.o serial_chimes_interface.o wrapper-C.o 
	ar rcs libchimes.a chimesFF.o serial_chimes_interface.o wrapper-C.o 

wrapper-F.o chimes.mod : wrapper-F.F90
	$(FCC) $(FFLAGS) -c wrapper-F.F90

clean:
	rm -f *.o *.x *.a *.mod
