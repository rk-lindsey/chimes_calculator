CC=g++
CCOMP=gcc
CCFLAGS = -g -std=c++11 
COMM1 = chimesFF.o serial_chimes_interface.o test_chimes_C_api.o wrapper-C.o

test_chimes_C_api : $(COMM1) 
	$(CC) $(CCFLAGS) $(COMM1) -o test_C.x

chimesFF.o : chimesFF.cpp
	$(CC) $(CCFLAGS) -c chimesFF.cpp

serial_chimes_interface.o : serial_chimes_interface.cpp
	$(CC) $(CCFLAGS) -c serial_chimes_interface.cpp

test_chimes_C_api.o : test_chimes_C_api.c
	$(CCOMP) -c test_chimes_C_api.c

wrapper-C.o : wrapper-C.cpp
	$(CC) -c wrapper-C.cpp

clean:
	rm -f *.o *.x
