FCC = gfortran
CC=g++
CCFLAGS = -fPIC -O3 -std=c++11 
FFLAGS = -O3 -fPIC
COMM1 = chimesFF.o serial_chimes_interface.o wrapper-F.o test_chimes_F_api.o wrapper-C.o 

test_F : $(COMM1) 
	$(FCC) $(FFLAGS) $(COMM1) -o test_F.x -lc++ 

chimesFF.o : chimesFF.cpp
	$(CC) $(CCFLAGS) -c chimesFF.cpp

serial_chimes_interface.o : serial_chimes_interface.cpp
	$(CC) $(CCFLAGS) -c serial_chimes_interface.cpp

test_chimes_F_api.o : test_chimes_F_api.F90
	$(FCC) $(FFLAGS) -c test_chimes_F_api.F90

wrapper-C.o : wrapper-C.cpp
	$(CC) $(CCFLAGS) -c wrapper-C.cpp

wrapper-F.o chimes.mod : wrapper-F.F90
	$(FCC) $(FFLAGS) -c wrapper-F.F90

clean:
	rm -f *.o *.x
